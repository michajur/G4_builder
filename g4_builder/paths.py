"""Construct run paths, deriving MD output names from the TPR filenames."""

from dataclasses import dataclass, field
from pathlib import Path
from string import Formatter
from types import MappingProxyType
import re

from .input_data import read_input_object
from .resources import data_path

_FILE_KEYS = {
    "strand", "full_ref", "reference", "plumed", "colvar", "topology", "posre",
    "conf", "box", "em_tpr", "tpr", "joined_trajectory", "frame", "log",
    "manifest",
    "em_mdp_template", "simulation_mdp_template",
}


def _format_filename(defaults, key, **values):
    """Expand a pattern, accepting only a filename without directories."""
    try:
        filename = defaults[key].format_map(values)
    except (KeyError, ValueError, AttributeError) as exc:
        raise ValueError(
            f"Invalid filename pattern for {key}: {defaults.get(key)!r}"
        ) from exc
    if (
        not filename
        or filename in {".", ".."}
        or Path(filename).name != filename
        or "\\" in filename
        or any(c.isspace() for c in filename)
    ):
        raise ValueError(
            f"{key} must resolve to a filename without directories "
            "or whitespace"
        )
    return filename


def load_file_names():
    """Load and validate package-level defaults; no per-run overrides."""
    defaults = read_input_object(data_path("file_names.json"))
    expected = _FILE_KEYS
    if set(defaults) != expected:
        raise ValueError(
            f"Invalid file_names.json keys; "
            f"missing={sorted(expected - set(defaults))}, "
            f"unknown={sorted(set(defaults) - expected)}"
        )
    for key in _FILE_KEYS:
        pattern = defaults[key]
        if not isinstance(pattern, str):
            raise ValueError(f"{key} must be a filename string")
        if key == "frame":
            allowed = {"index"}
            required = "index"
        elif key == "reference":
            allowed = {"name"}
            required = "name"
        else:
            allowed = set()
            required = None
        try:
            placeholders = {
                name
                for _, name, _, _ in Formatter().parse(pattern)
                if name is not None
            }
        except ValueError as exc:
            raise ValueError(f"Invalid filename pattern for {key}") from exc
        if placeholders - allowed:
            raise ValueError(
                f"Unsupported placeholders in {key}: "
                f"{sorted(placeholders - allowed)}"
            )
        if required and required not in placeholders:
            raise ValueError(
                f"{key} must include {{{required}}} to produce distinct files"
            )
        _format_filename(defaults, key, index=1, name="ref_1")
    for key in ("em_tpr", "tpr"):
        filename = Path(defaults[key])
        if filename.suffix != ".tpr" or not re.fullmatch(
            r"[A-Za-z0-9_-][A-Za-z0-9_.-]*", filename.stem
        ):
            raise ValueError(
                f"{key} must be a .tpr filename with a valid stem"
            )
    if defaults["em_tpr"] == defaults["tpr"]:
        raise ValueError(
            "Energy-minimization and simulation TPR filenames must differ"
        )
    return MappingProxyType(defaults)


def default_filename(key, **values):
    """Obtain a standalone component's default filename from the data file."""
    return _format_filename(load_file_names(), key, **values)


@dataclass(frozen=True)
class RunPaths:
    """Run directory and a snapshot of the package's filename defaults."""

    working_dir: Path
    _names: object = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        object.__setattr__(
            self, "working_dir", Path(self.working_dir).expanduser().resolve()
        )
        object.__setattr__(self, "_names", load_file_names())

    def path(self, name):
        """Place a filename or pattern inside the working directory."""
        return self.working_dir / name

    def filename(self, key, **values):
        """Expand a filename from this run's saved defaults."""
        return _format_filename(self._names, key, **values)

    def reference(self, name):
        """Return the output path for a named partial reference."""
        return self.path(self.filename("reference", name=name))

    def frame(self, index):
        """Return the output path for a selected frame's one-based rank."""
        return self.path(self.filename("frame", index=index))

    def to_dict(self):
        """Return absolute pathlib.Path values without scanning for files.

        Reference and frame paths retain their {name} and {index} placeholders.
        Prefixes are not paths and are omitted. Paths are included even if the
        output files do not exist.
        """
        names = (
            "working_dir", "strand", "full_ref", "plumed", "colvar",
            "topology",
            "posre", "conf", "box", "em_tpr", "em_structure", "mdp", "tpr",
            "simulation_structure", "trajectory", "joined_trajectory", "log",
            "manifest", "em_mdp_template", "simulation_mdp_template",
        )
        paths = {name: getattr(self, name) for name in names}
        paths["reference"] = self.path(self._names["reference"])
        paths["frame"] = self.path(self._names["frame"])
        return paths

    @property
    def em_prefix(self):
        return self.em_tpr.stem

    @property
    def simulation_prefix(self):
        return self.tpr.stem

    @property
    def strand(self):
        return self.path(self.filename("strand"))

    @property
    def full_ref(self):
        return self.path(self.filename("full_ref"))

    @property
    def plumed(self):
        return self.path(self.filename("plumed"))

    @property
    def colvar(self):
        return self.path(self.filename("colvar"))

    @property
    def topology(self):
        return self.path(self.filename("topology"))

    @property
    def posre(self):
        return self.path(self.filename("posre"))

    @property
    def conf(self):
        return self.path(self.filename("conf"))

    @property
    def box(self):
        return self.path(self.filename("box"))

    @property
    def em_tpr(self):
        return self.path(self.filename("em_tpr"))

    @property
    def em_structure(self):
        return self.em_tpr.with_suffix(".pdb")

    @property
    def mdp(self):
        return self.tpr.with_suffix(".mdp")

    @property
    def tpr(self):
        return self.path(self.filename("tpr"))

    @property
    def simulation_structure(self):
        return self.tpr.with_suffix(".pdb")

    @property
    def trajectory(self):
        return self.tpr.with_suffix(".xtc")

    @property
    def joined_trajectory(self):
        return self.path(self.filename("joined_trajectory"))

    @property
    def log(self):
        return self.path(self.filename("log"))

    @property
    def manifest(self):
        return self.path(self.filename("manifest"))

    @property
    def em_mdp_template(self):
        return data_path("mdp", self.filename("em_mdp_template"))

    @property
    def simulation_mdp_template(self):
        return data_path("mdp", self.filename("simulation_mdp_template"))


def require_files(*paths):
    """Raise FileNotFoundError for the first missing required file."""
    for path in paths:
        if not Path(path).is_file():
            raise FileNotFoundError(
                f"Required pipeline input not found: {path}"
            )
