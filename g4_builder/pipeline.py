"""Shared orchestration for Python callers and the terminal client."""

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Optional

from .config import SimulationConfig
from .extraction import TrajectoryFrameExtractor
from .gromacs import GromacsPreparation, SimulationOutputs
from .paths import RunPaths, require_files
from .pdb_io import read_pdb, save_pdb
from .protocol import Protocol, build_protocol
from .references import Reference, ReferenceSet, build_references
from .runner import GromacsRunner
from .strand import DNAStrandBuilder


@dataclass
class RunResult:
    """Paths, progress and results shared by all stages of one pipeline run."""

    paths: RunPaths
    protocol_steps: int
    full_ref_label: str
    simulation_steps: int
    reference_paths: dict = field(default_factory=dict)
    completed_stages: list = field(default_factory=list)
    outputs: Optional[SimulationOutputs] = None
    best_frames: dict = field(default_factory=dict)
    frame_paths: dict = field(default_factory=dict)


class FoldingPipeline:
    """Explicit stages: prepare, minimize, prepare_simulation, run, extract.

    Construction validates config but creates no files and starts no process.
    prepare() needs no GROMACS. Each successful stage updates run.json.
    GROMACS stdout/stderr go only through the injected runner's logger.
    """

    def __init__(self, config, runner=None):
        if not isinstance(config, SimulationConfig):
            raise TypeError("config must be SimulationConfig")
        self.config = config
        if runner is not None:
            self.runner = runner
        else:
            self.runner = GromacsRunner(
                config.working_dir, config.gromacs_executable
            )
        self.paths = self.runner.paths
        if self.paths.working_dir != config.working_dir:
            raise ValueError(
                "Runner and configuration working directories must match"
            )
        self.gromacs = GromacsPreparation(self.runner)
        self.result = None
        self.references = None
        self.protocol = None
        self.strand = None
        self.extractor = None

    def _require_prepared(self):
        if self.result is None:
            raise ValueError(
                "Call prepare() or load the pipeline with "
                "from_directory() first"
            )

    def _save_references(self):
        """Save references and copy their paths into the current result."""
        self.references.save_full_ref(self.paths.full_ref)
        self.references.save_partial_refs(self.paths.working_dir)
        self.result.reference_paths = self.references.saved_paths.copy()

    def _build_protocol(self, references, atom_serials):
        """Use the same protocol settings before and after pdb2gmx."""
        return build_protocol(
            references, atom_serials,
            schedule=self.config.schedule,
            print_stride=self.config.saving_stride,
            print_file=self.paths.colvar.name,
        )

    def _record(self, stage):
        """Write the run state atomically after a successful stage."""
        if stage is not None and stage not in self.result.completed_stages:
            self.result.completed_stages.append(stage)
        outputs = self.result.outputs
        metadata = {
            "schema_version": 1,
            # The saved location lets restoration relocate paths after a move.
            "config": dict(
                self.config.to_dict(),
                working_dir=str(self.config.working_dir),
            ),
            "protocol_steps": self.result.protocol_steps,
            "full_ref_label": self.result.full_ref_label,
            "simulation_steps": self.result.simulation_steps,
            "completed_stages": self.result.completed_stages,
            "references": [
                {
                    "name": entry.name,
                    "category": entry.category,
                    "resids": list(entry.resids),
                    "label": entry.label,
                    "is_full_core": entry.is_full_core,
                }
                for entry in self.references.entries
            ],
            "outputs": None if outputs is None else {
                "trajectories": [str(path) for path in outputs.trajectories],
                "structure": str(outputs.structure),
                "colvar": str(outputs.colvar),
            },
            "best_frames": self.result.best_frames,
            "frame_paths": {
                str(index): str(path)
                for index, path in self.result.frame_paths.items()
            },
        }
        temporary = self.paths.manifest.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(metadata, indent=2) + "\n")
        temporary.replace(self.paths.manifest)

    def save_metadata(self):
        """Persist configuration and results without marking a new stage."""
        self._require_prepared()
        self._record(None)

    def prepare(self, overwrite=False):
        """Build and save strand, references and protocol; no GROMACS needed.

        Refuse to replace an existing run unless overwrite=True. That option
        is only for input regeneration: it does not remove old MD outputs.
        Prefer a fresh working directory for a different simulation.
        """
        if not overwrite and any(
            path.exists()
            for path in (
                self.paths.manifest, self.paths.strand, self.paths.plumed
            )
        ):
            raise FileExistsError(
                "Run inputs already exist; use from_directory(), a new "
                "directory, or explicit overwrite=True"
            )
        builder = DNAStrandBuilder()
        strand = builder.build(self.config.sequence)
        references = build_references(
            self.config.polarity_pattern, self.config.helicity,
            strand, self.config.GT_resids, self.config.max_bulge_length,
        )
        protocol = self._build_protocol(references, strand.serial)
        # Complete validation and geometry before the first output file.
        self.paths.working_dir.mkdir(parents=True, exist_ok=True)
        self.strand = strand
        self.references = references
        self.protocol = protocol
        self.result = RunResult(
            self.paths,
            protocol.final_steps,
            protocol.full_ref_label,
            protocol.final_steps + self.config.after_folding_steps,
        )
        save_pdb(strand, self.paths.strand)
        self._save_references()
        protocol.write(self.paths.plumed)
        self._record("prepared")
        return self.result

    def minimize(self, pdb2gmx_timeout=60, discovery_timeout=2):
        """Create topology, remap changed serials, add a box and minimize."""
        self._require_prepared()
        if "simulation_prepared" in self.result.completed_stages:
            raise ValueError(
                "Use a new run directory to redo topology/minimization after "
                "simulation preparation"
            )
        self.gromacs.pdb2gmx(
            self.config.force_field,
            timeout=pdb2gmx_timeout,
            discovery_timeout=discovery_timeout,
        )
        structure = read_pdb(self.paths.conf)
        if self.references.needs_serial_remap(structure):
            print(
                "Reference atom order or serial numbering changed after "
                "pdb2gmx; "
                "remapping reference atom serials."
            )
            self.references.remap_serials(structure)
            self._save_references()
        self.protocol = self._build_protocol(self.references, structure.serial)
        self.protocol.write(self.paths.plumed)
        self.gromacs.add_box()
        self.gromacs.prep_enmin_file()
        self.gromacs.run_enmin()
        self._record("minimized")
        return self.result

    def prepare_simulation(self):
        """Prepare the MDP and TPR for the protocol plus post-folding steps."""
        self._require_prepared()
        if "minimized" not in self.result.completed_stages:
            raise ValueError("Call minimize() before prepare_simulation()")
        if "ran" in self.result.completed_stages:
            raise ValueError(
                "Simulation already ran; resume it with run() or use a "
                "new directory"
            )
        self.gromacs.prep_sim_file(
            self.result.simulation_steps,
            self.config.saving_stride,
            self.config.dielectric_const,
            dt_ps=self.config.dt_ps,
        )
        self._record("simulation_prepared")
        return self.result

    def run(self, timeout=None):
        """Run or resume the simulation and record its output files."""
        self._require_prepared()
        if "simulation_prepared" not in self.result.completed_stages:
            raise ValueError("Call prepare_simulation() before run()")
        self.gromacs.run_sim(timeout=timeout)
        self.result.outputs = self.gromacs.outputs
        self._record("ran")
        return self.result

    def extract(
        self, trajectories=None, structure=None, colvar=None, n=None,
        timeout=None,
    ):
        """Extract from recorded outputs or files supplied for an external run.

        Specify trajectories for an externally launched run. Default structure
        is the prepared TPR; no part number or executable is guessed.
        """
        self._require_prepared()
        outputs = self.result.outputs
        if trajectories is None:
            if outputs is None:
                raise ValueError(
                    "Run the simulation first, or supply trajectories from an "
                    "external run"
                )
            trajectories = outputs.trajectories
        if structure is None:
            structure = outputs.structure if outputs else self.paths.tpr
        if colvar is None:
            colvar = outputs.colvar if outputs else self.paths.colvar
        extractor = TrajectoryFrameExtractor(self.runner)
        extractor.find_best_frames(
            rmsd_label=self.result.full_ref_label,
            protocol_steps=self.result.protocol_steps,
            colvar=colvar,
            n=self.config.n_frames if n is None else n,
            dt_ps=self.config.dt_ps,
        )
        extractor.extract_frames(trajectories, structure, timeout=timeout)
        self.extractor = extractor
        self.result.best_frames = extractor.best_frames.copy()
        self.result.frame_paths = extractor.saved_paths.copy()
        self._record("extracted")
        return self.result

    @classmethod
    def from_directory(cls, directory, runner=None):
        """Restore a saved run in a fresh Python session, without GROMACS."""
        directory = Path(directory).expanduser().resolve()
        metadata = json.loads(RunPaths(directory).manifest.read_text())
        if metadata.get("schema_version") != 1:
            raise ValueError("Unsupported run manifest version")
        config = dict(metadata["config"])
        # Permit moving the whole run directory between machines.
        old_root = Path(config.pop("working_dir"))
        pipeline = cls(
            SimulationConfig.from_dict(config, working_dir=directory), runner
        )
        paths = pipeline.paths
        require_files(paths.strand, paths.full_ref, paths.plumed)
        entries = []
        for entry in metadata["references"]:
            entries.append(Reference(
                name=entry["name"],
                category=entry["category"],
                resids=tuple(entry["resids"]),
                label=entry["label"],
                atoms=read_pdb(paths.reference(entry["name"])),
                is_full_core=entry["is_full_core"],
            ))
        entries = tuple(entries)
        pipeline.references = ReferenceSet(read_pdb(paths.full_ref), entries)
        pipeline.references.saved_paths = {
            entry.name: paths.reference(entry.name) for entry in entries
        }
        pipeline.strand = read_pdb(paths.strand)
        pipeline.protocol = Protocol(
            paths.plumed.read_text(),
            metadata["protocol_steps"],
            metadata["full_ref_label"],
        )

        def relocate(value):
            """Relocate internal paths and retain external paths."""
            path = Path(value)
            if path.is_relative_to(old_root):
                return directory / path.relative_to(old_root)
            return path

        outputs = metadata["outputs"]
        if outputs:
            outputs = SimulationOutputs(
                tuple(relocate(path) for path in outputs["trajectories"]),
                relocate(outputs["structure"]),
                relocate(outputs["colvar"]),
            )
        pipeline.result = RunResult(
            paths=paths,
            protocol_steps=metadata["protocol_steps"],
            full_ref_label=metadata["full_ref_label"],
            simulation_steps=metadata["simulation_steps"],
            reference_paths=pipeline.references.saved_paths.copy(),
            completed_stages=list(metadata["completed_stages"]),
            outputs=outputs,
            best_frames={
                int(time): value
                for time, value in metadata["best_frames"].items()
            },
            frame_paths={
                int(index): relocate(path)
                for index, path in metadata["frame_paths"].items()
            },
        )
        return pipeline
