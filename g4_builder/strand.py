"""Nucleotide alignment and strand assembly."""
from pathlib import Path
import numpy as np
import pandas as pd
from .pdb_io import PDB_FIELDS, read_pdb, save_pdb
from .resources import data_path

class DNAStrandBuilder:
    """Load all 12 templates and assemble sequences containing at least two bases.

    ``templates`` maps template names (e.g. ``A5``, ``T``, ``C3``) to
    DataFrames. ``build`` returns the assembled DataFrame and stores it as
    ``strand``. Template DataFrames are never modified during assembly.

    Three non-collinear phosphate atoms determine a least-squares rigid
    alignment (Kabsch algorithm). Small template rounding differences can
    therefore leave a small residual alignment error.
    """

    _FIELDS = PDB_FIELDS
    _XYZ = ["x", "y", "z"]

    def __init__(self, nucleotide_dir=None):
        if nucleotide_dir is None:
            nucleotide_dir = data_path("nucleotides")
        self.nucleotide_dir = Path(nucleotide_dir)
        self.templates = {
            base + suffix: self.read_pdb(self.nucleotide_dir / f"{base}{suffix}.pdb")
            for base in "ATGC"
            for suffix in ("5", "", "3")
        }
        self.strand = None

    read_pdb = staticmethod(read_pdb)

    @classmethod
    def _anchors(cls, frame, residue):
        atoms = frame.loc[frame["res_seq"] == residue].copy()
        atoms["atom_name"] = atoms["atom_name"].replace({"O1P": "OP1", "O2P": "OP2"})
        anchors = []
        for name in ("P", "OP1", "OP2"):
            match = atoms.loc[atoms["atom_name"] == name, cls._XYZ]
            if len(match) != 1:
                raise ValueError(f"Residue {residue} must contain exactly one {name} anchor")
            anchors.append(match.iloc[0].to_numpy(dtype=float))
        anchors = np.asarray(anchors)
        if not np.isfinite(anchors).all() or np.linalg.matrix_rank(anchors - anchors.mean(axis=0)) < 2:
            raise ValueError("Phosphate anchors must be finite and non-collinear")
        return anchors

    def build(self, sequence, chain_id="A"):
        """Build a sequence written 5' to 3', returning its atom DataFrame.

        Whitespace and lowercase letters are accepted. A single nucleotide
        needs a template capped at both ends, which this template set lacks.
        """
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        sequence = "".join(sequence.split()).upper()
        if len(sequence) < 2:
            raise ValueError("At least two bases are required by the terminal templates")
        if set(sequence) - set("ATGC"):
            raise ValueError("sequence may contain only A, T, G, and C")
        if not isinstance(chain_id, str) or len(chain_id) != 1 or not chain_id.isascii() or not chain_id.isprintable():
            raise ValueError("chain_id must be one printable ASCII character")
        fragments = []
        target = None
        for index, base in enumerate(sequence):
            suffix = "5" if index == 0 else "3" if index == len(sequence) - 1 else ""
            frame = self.templates[base + suffix].copy(deep=True)
            if target is not None:
                source = self._anchors(frame, 1)
                source_center, target_center = source.mean(axis=0), target.mean(axis=0)
                u, _, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
                correction = np.eye(3)
                correction[-1, -1] = np.linalg.det(u @ vt)
                rotation = u @ correction @ vt
                frame[self._XYZ] = (frame[self._XYZ].to_numpy() - source_center) @ rotation + target_center
            if index < len(sequence) - 1:
                target = self._anchors(frame, 2)
            residue = frame.loc[frame["res_seq"] == 1].copy()
            if residue.empty:
                raise ValueError(f"Template {base + suffix} has no residue 1")
            residue["res_seq"] = index + 1
            residue["chain_id"] = chain_id
            fragments.append(residue)
        strand = pd.concat(fragments, ignore_index=True)
        strand["serial"] = np.arange(1, len(strand) + 1)
        self.strand = strand
        return strand

    def save_pdb(self, path):
        """Save the most recently built strand as fixed-width PDB records."""
        if self.strand is None:
            raise ValueError("Call build() before save_pdb()")
        save_pdb(self.strand, path)
