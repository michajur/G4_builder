"""GROMACS preparation stages, using an injected shared runner."""
from dataclasses import dataclass
from pathlib import Path
from numbers import Integral
import re
import tempfile

from .config import integer, positive
from .paths import require_files


@dataclass(frozen=True)
class SimulationOutputs:
    trajectories: tuple
    structure: Path
    colvar: Path


class GromacsPreparation:
    """Each method executes one stage and records its output paths.

    No force-field discovery or GROMACS execution occurs in __init__.
    runner is GromacsRunner (or an equivalent test double). Default file paths
    come from runner.paths; a separate paths argument is unnecessary.
    """

    def __init__(self, runner):
        self.runner = runner
        self.paths = runner.paths
        self.force_fields = {}
        self.saved_paths = {}
        self.dt_ps = None
        self.outputs = None

    def discover_force_fields(self, timeout=2, strand=None):
        """Discover the menu using an optional strand PDB path."""
        strand = self.paths.strand if strand is None else (self.paths.working_dir / Path(strand).expanduser()).resolve()
        require_files(strand)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = self.runner.probe(["pdb2gmx", "-f", strand,
                "-o", root / self.paths.conf.name, "-p", root / self.paths.topology.name,
                "-i", root / self.paths.posre.name], timeout)
        menu = output.partition("Select the Force Field:")[2].split("Select the Water Model:")[0]
        self.force_fields = {int(n): label.strip()
                             for n, label in re.findall(r"^\s*(\d+):\s*(.+)$", menu, re.MULTILINE)}
        return self.force_fields

    def pdb2gmx(self, force_field, strand=None, timeout=60, discovery_timeout=2):
        """Select a force field and prepare topology from a strand PDB.

        strand defaults to paths.strand. Relative paths use the runner's working
        directory; absolute paths and ~ are supported. Menu discovery uses the
        same input when a force-field menu number is supplied.
        """
        strand = self.paths.strand if strand is None else (self.paths.working_dir / Path(strand).expanduser()).resolve()
        require_files(strand)
        arguments = ["pdb2gmx", "-f", strand, "-o", self.paths.conf,
                     "-p", self.paths.topology, "-i", self.paths.posre, "-water", "tip3p"]
        answer = ""
        if isinstance(force_field, Integral) and not isinstance(force_field, bool):
            integer("force_field menu number", force_field, 1)
            if not self.force_fields:
                self.discover_force_fields(timeout=discovery_timeout, strand=strand)
            if force_field not in self.force_fields:
                raise ValueError(f"Force-field menu number {force_field} is not available; choices: {self.force_fields}")
            answer = f"{force_field}\n"
        elif isinstance(force_field, str) and re.fullmatch(r"[\w.+-]+", force_field):
            arguments += ["-ff", force_field]
        else:
            raise ValueError("Select a short force-field name or a discovered menu number")
        self.runner.run(arguments, input_text=answer, timeout=timeout)
        require_files(self.paths.conf, self.paths.topology)
        self.saved_paths.update(conf=self.paths.conf, topology=self.paths.topology)

    def add_box(self):
        require_files(self.paths.conf)
        self.runner.run(["editconf", "-f", self.paths.conf, "-bt", "dodecahedron",
                         "-d", "1.4", "-c", "-o", self.paths.box])
        require_files(self.paths.box)
        self.saved_paths["box"] = self.paths.box

    def prep_enmin_file(self):
        mdp = self.paths.em_mdp_template
        require_files(self.paths.box, self.paths.topology, mdp)
        self.runner.run(["grompp", "-f", mdp, "-p", self.paths.topology,
                         "-c", self.paths.box, "-o", self.paths.em_tpr, "-maxwarn", "1"])
        require_files(self.paths.em_tpr)
        self.saved_paths["em_tpr"] = self.paths.em_tpr

    def run_enmin(self):
        require_files(self.paths.em_tpr)
        self.runner.run(["mdrun", "-deffnm", self.paths.em_prefix,
                         "-c", self.paths.em_structure])
        require_files(self.paths.em_structure)
        self.saved_paths["em_structure"] = self.paths.em_structure

    def prep_sim_file(self, n_steps, saving_stride=5000, dielectric_const=15, dt_ps=0.002):
        """Render the packaged MDP and produce the centrally named .tpr."""

        # check input values 
        n_steps = integer("n_steps", n_steps, 1)
        saving_stride = integer("saving_stride", saving_stride, 1)
        dielectric_const = positive("dielectric_const", dielectric_const)
        dt_ps = positive("dt_ps", dt_ps)

        # check whether required files exist 
        require_files(self.paths.topology, self.paths.em_structure)

        # Read, modify and save MD parameter (.mdp) file 
        text = self.paths.simulation_mdp_template.read_text()
        replacements = {"SIM_STEPS": n_steps, "SAVE_STRIDE": saving_stride,
                        "EPSILON_VAL": dielectric_const, "DT_PS": dt_ps}
        for token, value in replacements.items():
            if token not in text:
                raise ValueError(f"Missing MDP template token: {token}")
            text = text.replace(token, str(value))
        self.paths.mdp.write_text(text)

        # Run gromacs command
        self.runner.run(["grompp", "-f", self.paths.mdp, "-p", self.paths.topology,
                         "-c", self.paths.em_structure, "-o", self.paths.tpr,
                         "-maxwarn", "2"])

        # Check if output was produced
        require_files(self.paths.tpr)

        # save dt_ps
        self.dt_ps = dt_ps

        # Collect output paths
        self.saved_paths.update(mdp=self.paths.mdp, tpr=self.paths.tpr)

    def run_sim(self, timeout=None):
        require_files(self.paths.tpr, self.paths.plumed)
        command = ["mdrun", "-deffnm", self.paths.simulation_prefix, "-cpi",
                   "-v", "-plumed", self.paths.plumed,
                   "-c", self.paths.simulation_structure]
        self.runner.run(command, timeout=timeout)

        require_files(self.paths.trajectory, self.paths.colvar)

        # The TPR is authoritative for atom ordering and works as trjconv -s.
        self.outputs = SimulationOutputs((self.paths.trajectory,), self.paths.tpr, self.paths.colvar)
