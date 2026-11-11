"""Load CLI simulation settings; execution inputs belong to command flags."""

from dataclasses import fields

from .config import SimulationConfig
from .input_data import read_input_object

# The Python API can change these fields. The CLI uses their defaults.
CLI_FIXED_CONFIG = {
    config_field.name: config_field.default
    for config_field in fields(SimulationConfig)
    if config_field.name == "dt_ps"
}


def validate_cli_config(config):
    """Require CLI defaults for settings configurable only through the API."""
    changed = [
        name for name, default in CLI_FIXED_CONFIG.items()
        if getattr(config, name) != default
    ]
    if changed:
        raise ValueError(
            "Run uses API-only settings differing from CLI defaults: "
            f"{', '.join(changed)}. Continue through the Python API or "
            "prepare a new run with CLI defaults."
        )


def load_input(path, working_dir):
    """Read simulation settings, taking working_dir from the CLI flag."""
    values = read_input_object(path)
    forbidden = set(values) & set(CLI_FIXED_CONFIG)
    if forbidden:
        raise ValueError(
            "These settings are fixed for CLI use and configurable only "
            f"through the Python API: {', '.join(sorted(forbidden))}"
        )
    return SimulationConfig.from_dict(values, working_dir=working_dir)
