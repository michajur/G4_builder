"""Validated simulation settings shared by the Python API and CLI."""

from dataclasses import dataclass, field, fields, MISSING
from collections.abc import Mapping
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Union
import math
import re

from .input_data import read_input_object


def integer(name, value, minimum=0):
    """Validate an integer setting, excluding bool, and return a Python int."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or value < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def positive(name, value):
    """Validate a positive finite real setting and return a Python float."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be positive and finite")
    return float(value)


@dataclass(frozen=True)
class ScheduleConfig:
    """Derive fixed restraint stages from four configurable MD step counts.

    GT hold/release stages follow the start of tetrad restraints.
    Core stages and tetrad stages use relative offsets.
    """
    gt_initial_steps: int = 400_000
    core_start: int = 1_000_000
    core_steps: int = 3_000_000
    tetrads_steps: int = 3_000_000

    def __post_init__(self):
        for name in (
            "gt_initial_steps", "core_start", "core_steps", "tetrads_steps"
        ):
            object.__setattr__(
                self, name, integer(name, getattr(self, name), 1)
            )
        if self.core_steps % 3:
            raise ValueError(
                "core_steps must be divisible by 3 so restraint stages use "
                "integer MD steps"
            )
        if self.core_start <= self.gt_initial_steps:
            raise ValueError(
                "core_start must be greater than gt_initial_steps"
            )
        if self.tetrads_steps <= self.core_steps // 3:
            raise ValueError(
                "tetrads_steps must be greater than core_steps / 3"
            )

    @property
    def gt_initial(self):
        """Initial G-tract ramp, measured from the simulation start."""
        return ((0, 0.0, 0), (self.gt_initial_steps, 0.0, 12000))

    @property
    def core(self):
        """Core restraint stages as offsets from each core block's start."""
        third = self.core_steps // 3
        return (
            (0, 0.0, 0),
            (third, 0.0, 20000),
            (2 * third, 0.0, 20000),
            (self.core_steps, 0.0, 0),
        )

    @property
    def core_spacing(self):
        """Core start interval, leaving one third of each block to overlap."""
        return 2 * (self.core_steps // 3)

    @property
    def tetrads(self):
        """Tetrad stages, with ramp duration tied to the core duration."""
        return (
            (0, 0.0, 0),
            (self.core_steps // 3, 0.0, 20000),
            (self.tetrads_steps, 0.0, 10000),
        )

    @property
    def gt_hold_kappa(self):
        """Continue the force constant reached by the last GT initial stage."""
        return self.gt_initial[-1][2]

    @property
    def gt_release_offset(self):
        """Reach zero GT kappa at the second tetrad stage's step offset."""
        return self.tetrads[1][0]


@dataclass(frozen=True)
class SimulationConfig:
    """Validated, immutable inputs shared by pipeline and CLI preparation."""

    sequence: str
    GT_resids: dict
    polarity_pattern: str
    helicity: str
    working_dir: Path
    force_field: Union[str, int]
    gromacs_executable: str = "gmx"
    dielectric_const: float = 15
    saving_stride: int = 5000
    after_folding_steps: int = 2_500_000
    dt_ps: float = 0.002
    max_bulge_length: int = 3
    n_frames: int = 5
    # Resolve the class when constructing the config, rather than capturing
    # its identity in a dataclass closure before a notebook module reload.
    schedule: ScheduleConfig = field(default_factory=lambda: ScheduleConfig())

    def __post_init__(self):
        # Validate sequence and G-tract correspondence before resolving paths.
        if not isinstance(self.sequence, str):
            raise ValueError("sequence must be a string")
        sequence = "".join(self.sequence.split()).upper()
        if len(sequence) < 2 or set(sequence) - set("ATGC"):
            raise ValueError(
                "sequence must contain at least two DNA bases (ATGC)"
            )
        if (
            not isinstance(self.polarity_pattern, str)
            or not self.polarity_pattern
            or set(self.polarity_pattern) - set("as")
        ):
            raise ValueError("polarity_pattern must contain only 'a' and 's'")
        if self.helicity not in ("RH", "LH"):
            raise ValueError("helicity must be 'RH' or 'LH'")
        if not isinstance(self.GT_resids, Mapping) or len(self.GT_resids) != 4:
            raise ValueError("Exactly four G-tracts are required for tetrads")
        tracts = {}
        seen = set()
        for name, values in self.GT_resids.items():
            if not isinstance(name, str) or not name:
                raise ValueError("G-tract names must be non-empty strings")
            if isinstance(values, (set, frozenset, dict, str, bytes)):
                raise ValueError(
                    "Each G-tract must be an ordered residue sequence"
                )
            try:
                residues = tuple(
                    integer("residue number", value, 1) for value in values
                )
            except TypeError as exc:
                raise ValueError(
                    "Each G-tract must be an ordered residue sequence"
                ) from exc
            if (
                len(residues) < 2
                or len(residues) != len(self.polarity_pattern)
            ):
                raise ValueError(
                    "Each G-tract must have at least two residues and match "
                    "polarity_pattern length"
                )
            for number in residues:
                if number > len(sequence) or sequence[number - 1] != "G":
                    raise ValueError(
                        f"Residue {number} must exist and be guanine"
                    )
                if number in seen:
                    raise ValueError(f"Duplicate G-tract residue: {number}")
                seen.add(number)
            tracts[name] = residues

        # Force fields accept a GROMACS name or an interactive menu index.
        if isinstance(self.force_field, Integral) and not isinstance(
            self.force_field, bool
        ):
            object.__setattr__(
                self, "force_field",
                integer("force_field menu number", self.force_field, 1),
            )
        elif not isinstance(self.force_field, str) or not re.fullmatch(
            r"[\w.+-]+", self.force_field
        ):
            raise ValueError(
                "force_field must be a short -ff name or a positive integer "
                "menu number"
            )
        if not str(self.gromacs_executable).strip():
            raise ValueError("gromacs_executable must not be empty")
        for name in ("saving_stride", "n_frames"):
            integer(name, getattr(self, name), 1)
        for name in ("after_folding_steps", "max_bulge_length"):
            integer(name, getattr(self, name))
        positive("dielectric_const", self.dielectric_const)
        positive("dt_ps", self.dt_ps)
        if not isinstance(self.schedule, ScheduleConfig):
            raise ValueError("schedule must be ScheduleConfig")

        # Store the remaining normalized inputs after all validation succeeds.
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "GT_resids", MappingProxyType(tracts))
        object.__setattr__(
            self, "working_dir", Path(self.working_dir).expanduser().resolve()
        )
        executable = str(self.gromacs_executable)
        if "/" in executable:
            executable = str(Path(executable).expanduser().resolve())
        object.__setattr__(self, "gromacs_executable", executable)

    def to_dict(self):
        """Return JSON-compatible input settings, excluding working_dir."""
        from dataclasses import asdict, fields

        result = {}
        for config_field in fields(self):
            if config_field.name != "working_dir":
                result[config_field.name] = getattr(self, config_field.name)
        result["GT_resids"] = {
            key: list(values) for key, values in self.GT_resids.items()
        }
        result["schedule"] = asdict(self.schedule)
        return result

    @classmethod
    def from_json(cls, path, working_dir):
        """Load simulation settings from JSON and a separate working directory.

        Use the same parameter names as the Python constructor. API-only
        settings such as dt_ps are allowed. Omitted options use their defaults.
        working_dir is required as an argument and must not occur in the JSON.
        Relative paths use the caller's cwd.

        Only SimulationConfig fields are accepted. File paths and execution
        controls for individual stages are passed to those methods separately.
        """
        return cls.from_dict(read_input_object(path), working_dir=working_dir)

    @classmethod
    def from_dict(cls, values, working_dir):
        """Load settings with working_dir supplied only as a separate argument.

        Reject working_dir in values. Copy the input mapping before validation;
        missing optional settings use their defaults.
        """
        if not isinstance(values, Mapping):
            raise ValueError("Simulation settings must be a mapping")
        if "working_dir" in values:
            raise ValueError(
                "working_dir must be supplied separately, not inside the "
                "input settings"
            )
        if working_dir is None:
            raise ValueError("working_dir must be provided separately")
        values = dict(values)
        values["working_dir"] = working_dir
        config_fields = {f.name: f for f in fields(cls)}
        unknown = set(values) - set(config_fields)
        if unknown:
            raise ValueError(f"Unknown input settings: {sorted(unknown)}")
        required = set()
        for name, config_field in config_fields.items():
            if (
                config_field.default is MISSING
                and config_field.default_factory is MISSING
            ):
                required.add(name)
        missing = required - set(values)
        if missing:
            raise ValueError(
                f"Missing mandatory input settings: {sorted(missing)}"
            )
        try:
            if isinstance(values.get("schedule"), dict):
                values["schedule"] = ScheduleConfig(**values["schedule"])
            return cls(**values)
        except (TypeError, AttributeError) as exc:
            raise ValueError(f"Invalid input configuration: {exc}") from exc
