"""Fixed-width PDB atom tables, independent of strand construction."""

from pathlib import Path

import numpy as np
import pandas as pd

PDB_FIELDS = (
    ("record_name", 0, 6, str),
    ("serial", 6, 11, int),
    ("atom_name", 12, 16, str),
    ("alt_loc", 16, 17, str),
    ("res_name", 17, 20, str),
    ("chain_id", 21, 22, str),
    ("res_seq", 22, 26, int),
    ("i_code", 26, 27, str),
    ("x", 30, 38, float),
    ("y", 38, 46, float),
    ("z", 46, 54, float),
    ("occupancy", 54, 60, float),
    ("temp_factor", 60, 66, float),
    ("seg_id", 72, 76, str),
    ("element", 76, 78, str),
    ("charge", 78, 80, str),
)


def read_pdb(path):
    """Read fixed-width ATOM/HETATM records; ignore non-atom records.

    Optional blank numeric fields are represented by NaN. Atom names
    are stripped; their standard PDB alignment is restored on export.
    """
    rows = []
    with Path(path).expanduser().open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            row = {}
            try:
                for name, start, end, conversion in PDB_FIELDS:
                    value = line[start:end].strip()
                    if conversion is float and not value and start >= 54:
                        row[name] = np.nan
                    else:
                        row[name] = conversion(value)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid PDB field in {path}, line {line_number}"
                ) from exc
            rows.append(row)
    if not rows:
        raise ValueError(f"No atom records found in {path}")
    return pd.DataFrame(rows, columns=[field[0] for field in PDB_FIELDS])


def save_pdb(df, path, add_ter_end=True):
    """Save an atom DataFrame as fixed-width PDB records.

    Expect the columns produced by read_pdb(). Preserve row order and atom
    and residue numbering. Set add_ter_end=False to omit TER and END records.
    Accept a string or Path for the output path. The DataFrame is unchanged.
    """
    lines = []
    for atom in df.itertuples(index=False):
        name = atom.atom_name
        # One-letter elements start in column 14, except digit-led names.
        if len(name) < 4 and not name[0].isdigit() and len(atom.element) == 1:
            name = " " + name
        if pd.isna(atom.occupancy):
            occupancy = " " * 6
        else:
            occupancy = f"{atom.occupancy:6.2f}"
        if pd.isna(atom.temp_factor):
            temperature = " " * 6
        else:
            temperature = f"{atom.temp_factor:6.2f}"
        line = (
            f"{atom.record_name:<6}{atom.serial:5d} {name:<4}{atom.alt_loc:1}"
            f"{atom.res_name:>3} {atom.chain_id:1}"
            f"{atom.res_seq:4d}{atom.i_code:1}   "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{occupancy}{temperature}"
            f"      {atom.seg_id:<4}{atom.element:>2}{atom.charge:>2}"
        )
        if len(line) != 80:
            raise ValueError(
                "DataFrame contains values exceeding fixed-width PDB fields"
            )
        lines.append(line)
    if add_ter_end:
        lines.extend(["TER", "END"])
    Path(path).expanduser().write_text("".join(line + "\n" for line in lines))
