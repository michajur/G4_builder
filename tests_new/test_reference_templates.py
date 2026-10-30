"""Geometry, template-source selection and reference integration."""
from pathlib import Path
import json
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from g4_builder import DNAStrandBuilder, build_reference_template, build_references
from g4_builder.pdb_io import read_pdb
from g4_builder.resources import data_path
from tests_new.test_library import config


class ReferenceTemplateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_geometry_and_residue_order_with_distinct_pair_properties(self):
        layers = data_path("reference_templates", "layers")
        original_bytes = {base: (layers / f"{base}.pdb").read_bytes() for base in "as"}
        templates = {base: read_pdb(layers / f"{base}.pdb") for base in "as"}
        properties = {pair: {"rise": 2.1 + i, "twist_RH": 11 + 7 * i,
                             "twist_LH": -13 - 9 * i}
                      for i, pair in enumerate(("aa", "as", "sa", "ss"))}
        properties_path = self.root / "stack_properties.json"
        properties_path.write_text(json.dumps(properties))
        def resource(*parts):
            return properties_path if parts == ("stack_properties.json",) else data_path(*parts)
        with patch("g4_builder.template_builder.data_path", side_effect=resource):
            for pattern in ("a", "saa", "as", "sass"):
                for helicity in ("RH", "LH"):
                    with self.subTest(pattern=pattern, helicity=helicity):
                        template = build_reference_template(pattern, helicity)
                        self.assertEqual(len(template), 44 * len(pattern))
                        self.assertEqual(template.serial.tolist(), list(range(1, len(template) + 1)))
                        self.assertEqual(template.res_seq.drop_duplicates().tolist(),
                                         list(range(1, 4 * len(pattern) + 1)))
                        # Independently compose parent-relative 4x4 transforms.
                        transform = np.eye(4)
                        for layer_index, base in enumerate(pattern):
                            if layer_index:
                                pair = properties[pattern[layer_index - 1:layer_index + 1]]
                                angle = np.deg2rad(pair[f"twist_{helicity}"])
                                c, s = np.cos(angle), np.sin(angle)
                                transform = transform @ np.array([
                                    [c, -s, 0, 0], [s, c, 0, 0],
                                    [0, 0, 1, pair["rise"]], [0, 0, 0, 1]])
                            for residue in range(1, 5):
                                atoms = templates[base].loc[templates[base].res_seq.eq(residue)]
                                number = (residue - 1) * len(pattern) + layer_index + 1
                                actual = template.loc[template.res_seq.eq(number)]
                                self.assertEqual(actual.atom_name.tolist(), atoms.atom_name.tolist())
                                xyz = atoms[["x", "y", "z"]].to_numpy()
                                expected = (transform @ np.column_stack((xyz, np.ones(len(xyz)))).T).T[:, :3]
                                np.testing.assert_allclose(actual[["x", "y", "z"]], expected, atol=1e-12)
        for base in "as":
            self.assertEqual((layers / f"{base}.pdb").read_bytes(), original_bytes[base])

    def test_loading_remains_default_and_building_supports_missing_stored_template(self):
        c = config(self.root)
        strand = DNAStrandBuilder().build(c.sequence)
        original = strand.copy(deep=True)
        with patch("g4_builder.references.build_reference_template",
                   side_effect=AssertionError("Loading must not invoke the builder")):
            loaded = build_references(c.polarity_pattern, c.helicity, strand, c.GT_resids)
        with self.assertRaises(FileNotFoundError):
            build_references(c.polarity_pattern, "LH", strand, c.GT_resids)

        built = build_references(c.polarity_pattern, "LH", strand, c.GT_resids, build_template=True)
        template = build_reference_template(c.polarity_pattern, "LH")
        for tract_index, resids in enumerate(c.GT_resids.values()):
            for layer_index, resid in enumerate(resids):
                number = tract_index * len(c.polarity_pattern) + layer_index + 1
                expected = template.loc[template.res_seq.eq(number)]
                actual = built.full_ref.loc[built.full_ref.res_seq.eq(resid)]
                self.assertEqual(actual.atom_name.tolist(), expected.atom_name.tolist())
                np.testing.assert_allclose(actual[["x", "y", "z"]], expected[["x", "y", "z"]])
        self.assertEqual(built.ref_info, loaded.ref_info)
        self.assertEqual(built.full_ref_label, loaded.full_ref_label)
        self.assertTrue(built.full_ref.temp_factor.eq(1.0).all())
        self.assertTrue(built.full_ref.chain_id.eq("").all())
        self.assertTrue(built.full_ref.element.eq("").all())
        pd.testing.assert_frame_equal(strand, original)

    def test_custom_template_root_uses_its_layers_without_stored_template(self):
        layers = self.root / "layers"
        shutil.copytree(data_path("reference_templates", "layers"), layers)
        # Distinguish this custom source from the packaged layer coordinates.
        path = layers / "a.pdb"
        lines = path.read_text().splitlines(keepends=True)
        lines[0] = lines[0][:30] + f"{42.0:8.3f}" + lines[0][38:]
        path.write_text("".join(lines))
        expected = build_reference_template("saa", "RH", layer_dir=layers)
        self.assertEqual(expected.loc[expected.res_seq.eq(2), "x"].iloc[0],
                         42.0 * np.cos(np.deg2rad(25))
                         - read_pdb(path).y.iloc[0] * np.sin(np.deg2rad(25)))
        c = config(self.root)
        strand = DNAStrandBuilder().build(c.sequence)
        with self.assertRaises(FileNotFoundError):
            build_references("saa", "RH", strand, c.GT_resids, template_dir=self.root)
        refs = build_references("saa", "RH", strand, c.GT_resids,
                                template_dir=self.root, build_template=True)
        actual = refs.full_ref.loc[refs.full_ref.res_seq.eq(c.GT_resids["GT_1"][1])]
        np.testing.assert_allclose(actual[["x", "y", "z"]],
                                   expected.loc[expected.res_seq.eq(2), ["x", "y", "z"]])
        self.assertFalse((self.root / "RH_3tet").exists())


if __name__ == "__main__":
    unittest.main()
