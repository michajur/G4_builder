"""Build G-quadruplex reference templates from flat a/s layers."""

from pathlib import Path

import numpy as np
import pandas as pd

from .pdb_io import read_pdb
from .input_data import read_input_object
from .resources import data_path


def build_reference_template(polarity_pattern, helicity, layer_dir=None):
    """Stack flat a/s layers from bottom to top and return a PDB DataFrame.

    Read a.pdb and s.pdb from packaged reference_templates/layers, or from
    layer_dir if supplied. Each layer must contain residues 1, 2, 3, 4
    with finite coordinates in the Z=0 plane.

    For each adjacent pair, use its rise (angstroms) and twist_RH or twist_LH
    (degrees) from data/stack_properties.json, read on each call. Twists are
    signed rotations about the Z axis through the XY origin; positive angles
    rotate counterclockwise when viewed from positive Z. Rise and twist
    accumulate relative to the preceding layer. The bottom layer retains its
    original coordinates at Z=0.

    Group residues by their original number, ordered bottom to top within each
    group. For three layers, residues 1-3 come from original residue 1,
    residues 4-6 from original residue 2, and so on. Keep atom order within
    residues and renumber atoms consecutively. No files are modified or saved.
    A single-letter pattern returns one layer.
    """
    if (
        not isinstance(polarity_pattern, str)
        or not polarity_pattern
        or set(polarity_pattern) - {"a", "s"}
    ):
        raise ValueError(
            "polarity_pattern must be a nonempty string containing "
            "only 'a' and 's'"
        )
    if helicity not in ("RH", "LH"):
        raise ValueError("helicity must be 'RH' or 'LH'")

    stack_properties = read_input_object(data_path("stack_properties.json"))
    if layer_dir is None:
        root = data_path("reference_templates", "layers")
    else:
        root = Path(layer_dir).expanduser().resolve()
    templates = {base: read_pdb(root / f"{base}.pdb") for base in "as"}
    for base, template in templates.items():
        xyz = template[["x", "y", "z"]].to_numpy(dtype=float)
        if (
            not np.isfinite(xyz).all()
            or not np.allclose(xyz[:, 2], 0, atol=1e-8)
        ):
            raise ValueError(
                f"{base}.pdb must contain finite coordinates in the Z=0 plane"
            )
        if set(template["res_seq"]) != {1, 2, 3, 4}:
            raise ValueError(
                f"{base}.pdb must contain exactly residues 1, 2, 3, 4"
            )

    layers = []
    height = 0.0
    angle_degrees = 0.0
    twist_key = f"twist_{helicity}"
    for index, base in enumerate(polarity_pattern):
        if index > 0:
            # Each transform is relative to the previous layer, not the bottom.
            pair = polarity_pattern[index - 1:index + 1]
            settings = stack_properties[pair]
            rise = float(settings["rise"])
            twist = float(settings[twist_key])
            if not np.isfinite([rise, twist]).all():
                raise ValueError(
                    f"Rise and {twist_key} for {pair!r} must be finite"
                )
            height += rise
            angle_degrees += twist

        layer = templates[base].copy(deep=True)
        angle = np.deg2rad(angle_degrees)
        cosine, sine = np.cos(angle), np.sin(angle)
        # Row-vector coordinates require the transpose of the usual Rz matrix.
        rotation = np.array([[cosine, sine], [-sine, cosine]])
        layer[["x", "y"]] = layer[["x", "y"]].to_numpy() @ rotation
        layer["z"] = layer["z"] + height

        # Number residues consecutively along each column through the stack.
        layer["res_seq"] = (
            (layer["res_seq"] - 1) * len(polarity_pattern) + index + 1
        )
        layers.append(layer)

    stack = pd.concat(layers, ignore_index=True)
    stack = stack.sort_values("res_seq", kind="stable").reset_index(drop=True)
    stack["serial"] = np.arange(1, len(stack) + 1)
    return stack
