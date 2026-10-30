"""Generate PLUMED text without writing files or launching commands."""
from dataclasses import dataclass
from pathlib import Path
from .config import ScheduleConfig, integer
from .paths import default_filename


@dataclass(frozen=True)
class Protocol:
    text: str
    final_steps: int
    full_ref_label: str

    def write(self, path):
        Path(path).expanduser().write_text(self.text)


def build_protocol(references, atom_serials, schedule=None, print_stride=5000,
                   print_file=None):
    """Order GT/core/tetrad restraints, referencing PDBs by filename.

    Reference PDBs, the PLUMED file and simulation TPR must share the GROMACS
    working directory when running the simulation.
    """
    schedule = schedule or ScheduleConfig()
    print_file = default_filename("colvar") if print_file is None else print_file
    integer("print_stride", print_stride, 1)
    if not str(print_file) or any(c.isspace() for c in str(print_file)):
        raise ValueError("PLUMED path tokens must be non-empty and contain no whitespace")
    serials = sorted(integer("atom serial", n, 1) for n in atom_serials)
    if not serials or len(serials) != len(set(serials)):
        raise ValueError("Atom serials must be non-empty and unique")
    entity = (f"{serials[0]}-{serials[-1]}" if serials[-1] - serials[0] + 1 == len(serials)
              else ",".join(map(str, serials)))
    lines = ["RESTART", "", f"WHOLEMOLECULES ENTITY0={entity}", ""]
    for entry in references.entries:
        lines.append(f"RMSD REFERENCE={entry.filename} TYPE=OPTIMAL LABEL={entry.label}")
    lines.append("")
    final_steps = 0

    def restraint(labels, stages, start=0):
        nonlocal final_steps
        if not labels:
            return
        repeated = lambda value: ",".join([str(value)] * len(labels))
        lines.extend(["MOVINGRESTRAINT ...", f"  ARG={','.join(labels)}", f"  VERSE={repeated('U')}"])
        for i, (offset, target, kappa) in enumerate(stages):
            step = start + offset
            final_steps = max(final_steps, step)
            lines.append(f"  STEP{i}={step}  AT{i}={repeated(target)}  KAPPA{i}={repeated(kappa)}")
        lines.extend(["... MOVINGRESTRAINT", ""])

    groups = {name: [e.label for e in references.entries if e.category == name]
              for name in ("GT", "core", "tetrads")}
    tetrad_start = schedule.core_start + len(groups["core"]) * schedule.core_spacing
    gt_stages = schedule.gt_initial + ((tetrad_start, 0.0, schedule.gt_hold_kappa),
                                       (tetrad_start + schedule.gt_release_offset, 0.0, 0))
    restraint(groups["GT"], gt_stages)
    for index, label in enumerate(groups["core"]):
        restraint([label], schedule.core, schedule.core_start + index * schedule.core_spacing)
    restraint(groups["tetrads"], schedule.tetrads, tetrad_start)
    labels = ",".join(entry.label for entry in references.entries)
    lines.append(f"PRINT ARG={labels} STRIDE={print_stride} FILE={print_file}")
    return Protocol("\n".join(lines) + "\n", final_steps, references.full_ref_label)
