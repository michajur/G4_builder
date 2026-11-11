"""Locations of templates distributed with the package."""

from importlib.resources import files
from pathlib import Path


def data_path(*parts):
    """Locate packaged data in a source tree or installed wheel."""
    path = Path(str(files("g4_builder").joinpath("data", *parts)))
    if not path.exists():
        raise FileNotFoundError(f"Package template not found: {path}")
    return path
