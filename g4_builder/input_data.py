"""Shared JSON parsing for API and CLI callers."""
import json
from pathlib import Path

def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_input_object(path):
    with Path(path).expanduser().open() as handle:
        values = json.load(handle, object_pairs_hook=_unique_object)
    if not isinstance(values, dict):
        raise ValueError("Input JSON must contain an object")
    return values


