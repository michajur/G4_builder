"""Shared JSON parsing for API and CLI callers."""

import json
from pathlib import Path


def _unique_object(pairs):
    """Build a JSON object without silently accepting duplicate keys."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_input_object(path):
    """Load a JSON object, expanding the user's home directory in its path."""
    with Path(path).expanduser().open() as handle:
        values = json.load(handle, object_pairs_hook=_unique_object)
    if not isinstance(values, dict):
        raise ValueError("Input JSON must contain an object")
    return values
