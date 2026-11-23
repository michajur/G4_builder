"""Reference construction and explicit metadata for each restraint."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import integer
from .pdb_io import save_pdb
from .paths import default_filename
from .template_builder import build_reference_template

_G_BASE_ATOMS = (
    "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"
)
_XYZ = ["x", "y", "z"]


@dataclass(frozen=True)
class Reference:
    """One restraint's selected atoms, residue membership and RMSD label."""

    name: str
    category: str
    resids: tuple
    label: str
    atoms: pd.DataFrame
    is_full_core: bool = False

    @property
    def filename(self):
        """Filename used by both PDB saving and PLUMED generation."""
        return default_filename("reference", name=self.name)


@dataclass
class ReferenceSet:
    """Built reference atoms and metadata, with saving and serial remapping."""

    full_ref: pd.DataFrame
    entries: tuple
    saved_paths: dict = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    @property
    def partial_refs(self):
        """Map partial-reference names to their atom DataFrames."""
        return {entry.name: entry.atoms for entry in self.entries}

    @property
    def ref_info(self):
        """Group ordered residue lists by the restraint type they describe."""
        groups = {}
        for group in ("GT", "core", "tetrads"):
            groups[group] = []
            for entry in self.entries:
                if entry.category == group:
                    groups[group].append(list(entry.resids))
        return groups

    @property
    def full_ref_label(self):
        """RMSD label of the last cumulative core reference."""
        return next(
            entry.label for entry in self.entries if entry.is_full_core
        )

    def save_full_ref(self, path):
        """Save the full reference without TER or END records; return None."""
        save_pdb(self.full_ref, path, add_ter_end=False)

    def save_partial_refs(self, directory):
        """Save partial references and print/store their paths."""
        root = Path(directory).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)
        self.saved_paths = {}
        for entry in self.entries:
            path = root / entry.filename
            save_pdb(entry.atoms, path, add_ter_end=False)
            self.saved_paths[entry.name] = path
            print(f"Saved {entry.name}: {path}")

    def _serial_updates(self, structure):
        """Validate atom identities and collect serials without mutation."""
        keys = ["res_seq", "atom_name"]
        mapping = structure.set_index(keys)["serial"]
        if mapping.index.has_duplicates:
            raise ValueError(
                "Processed structure has ambiguous residue/atom identifiers"
            )
        tables = [self.full_ref, *(entry.atoms for entry in self.entries)]
        updates = []
        for table in tables:
            try:
                values = mapping.loc[
                    pd.MultiIndex.from_frame(table[keys])
                ].to_numpy()
            except KeyError as exc:
                raise ValueError(
                    "Processed structure is missing a reference base atom"
                ) from exc
            updates.append((table, values))
        return updates

    def needs_serial_remap(self, structure):
        """Check whether any reference atom's serial changed after pdb2gmx.

        Match by residue number and atom name, not by DataFrame row position.
        Check the full reference and every partial reference. Missing or
        ambiguous atoms raise ValueError even when no remapping is needed.
        """
        updates = self._serial_updates(structure)
        return any(
            not np.array_equal(table["serial"].to_numpy(), values)
            for table, values in updates
        )

    def remap_serials(self, structure):
        """Update reference serials by residue number and heavy-base atom name.

        Coordinates and restraint ordering are retained. All matches are
        validated before any table is modified.
        """
        for table, values in self._serial_updates(structure):
            table["serial"] = values


def build_references(
    polarity_pattern, helicity, strand, GT_resids,
    max_bulge_length=3, layer_dir=None,
):
    """Build a full reference and GT/core/tetrad partials as a ReferenceSet.

    G-tract dictionary order and residue order define template correspondence.
    Generate the template using build_reference_template() from the packaged
    flat a.pdb/s.pdb layers, or from layer_dir when supplied. Rise and twist
    values come from data/stack_properties.json. Atom names must match the
    template row for row before coordinates are copied. Input tables and
    G-tracts are unchanged; no files are written.
    Use the returned object's save_full_ref(), save_partial_refs(),
    needs_serial_remap() and remap_serials() methods for subsequent steps.
    """

    # remove it later 
    print(
        "!!! WARNING !!!\n"\
        "Reference templates are built based on ad-hoc parameters"\
        "(Twist and rise of tetrads) so gemotry of the G-core might\n"\
        "be a little off. Parameters will be derived from experimental"\
        "structures and updated in the future."
    )

    if (
        not isinstance(polarity_pattern, str)
        or not polarity_pattern
        or set(polarity_pattern) - set("as")
    ):
        raise ValueError("polarity_pattern must contain 'a' and 's'")
    if helicity not in ("RH", "LH"):
        raise ValueError("helicity must be RH or LH")
    max_bulge_length = integer("max_bulge_length", max_bulge_length)
    if not isinstance(strand, pd.DataFrame):
        raise TypeError("strand must be a DataFrame")
    tracts = {}
    for name, values in GT_resids.items():
        tracts[name] = tuple(
            integer("residue", number, 1) for number in values
        )
    if len(tracts) != 4 or any(
        len(values) != len(polarity_pattern) for values in tracts.values()
    ):
        raise ValueError(
            "Four G-tracts, each matching polarity_pattern length, "
            "are required"
        )
    resids = []
    for values in tracts.values():
        resids.extend(values)
    if len(set(resids)) != len(resids):
        raise ValueError("G-tract residues must be unique")
    template = build_reference_template(
        polarity_pattern, helicity, layer_dir=layer_dir
    )

    # Assign coordinates in G-tract order, then sort the full reference.
    fragments = []
    for values in tracts.values():
        for number in values:
            residue = strand.loc[strand.res_seq == number]
            if residue.empty or not residue.res_name.isin(["DG", "G"]).all():
                raise ValueError(f"Residue {number} must exist and be guanine")
            atoms = residue.loc[residue.atom_name.isin(_G_BASE_ATOMS)]
            if (
                len(atoms) != len(_G_BASE_ATOMS)
                or set(atoms.atom_name) != set(_G_BASE_ATOMS)
            ):
                raise ValueError(
                    f"Residue {number} must contain each guanine "
                    "base atom once"
                )
            fragments.append(atoms)
    reference = pd.concat(fragments, ignore_index=True).copy()
    if len(reference) != len(template):
        raise ValueError("Reference and template atom counts differ")
    if not np.array_equal(
        reference.atom_name.to_numpy(), template.atom_name.to_numpy()
    ):
        raise ValueError(
            "Reference and template atom names do not match row for row"
        )
    coordinates = template[_XYZ].to_numpy(dtype=float)
    if not np.isfinite(coordinates).all():
        raise ValueError("Reference template coordinates must be finite")
    reference[_XYZ] = coordinates
    reference["temp_factor"] = 1.0
    reference["chain_id"] = ""
    reference["element"] = ""
    full_ref = reference.sort_values(
        "res_seq", kind="stable"
    ).reset_index(drop=True)

    # Split G-tracts at long bulges and assemble the partial reference groups.
    blocks = []
    for values in tracts.values():
        block = [values[0]]
        for previous, current in zip(values, values[1:]):
            if abs(current - previous) <= max_bulge_length + 1:
                block.append(current)
            else:
                blocks.append(block)
                block = [current]
        blocks.append(block)
    blocks.sort(key=min)
    gt = [tuple(block) for block in blocks if len(block) > 1]
    core = []
    combined = list(blocks[0])
    for block in blocks[1:]:
        combined = combined + block
        core.append(tuple(combined))
    tetrads = list(zip(*tracts.values()))
    entries = []
    for group, residues in (("GT", gt), ("core", core), ("tetrads", tetrads)):
        for position, resids in enumerate(residues):
            index = len(entries) + 1
            atoms = full_ref.loc[
                full_ref.res_seq.isin(resids)
            ].copy().reset_index(drop=True)
            # The last cumulative core contains every block by construction.
            is_full_core = group == "core" and position == len(core) - 1
            entries.append(Reference(
                f"ref_{index}", group, tuple(resids), f"rm{index}",
                atoms, is_full_core,
            ))
    return ReferenceSet(full_ref, tuple(entries))
