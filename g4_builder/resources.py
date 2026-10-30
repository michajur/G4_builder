"""Locations of templates distributed with the package."""
from importlib.resources import files
from pathlib import Path


def data_path(*parts):
    """Return a resource path from a normal source or installed-wheel layout."""
    path = Path(str(files("g4_builder").joinpath("data", *parts)))
    if not path.exists():
        raise FileNotFoundError(f"Package template not found: {path}")
    return path
