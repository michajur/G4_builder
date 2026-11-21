"""Generated template geometry, custom layers and reference integration."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import json
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from g4_builder import (
    DNAStrandBuilder, FoldingPipeline, build_reference_template, build_references,
    configure_logging,
)
from g4_builder.cli import main
from g4_builder.pdb_io import read_pdb
from g4_builder.resources import data_path
from tests_new.test_library import config


class ReferenceTemplateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.addCleanup(lambda: configure_logging(console=False))

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

    def test_references_always_generate_templates_from_flat_layers(self):
        c = config(self.root)
        strand = DNAStrandBuilder().build(c.sequence)
        original = strand.copy(deep=True)
        with patch("g4_builder.references.build_reference_template",
                   wraps=build_reference_template) as builder, patch(
                       "g4_builder.template_builder.read_pdb", wraps=read_pdb
                   ) as reader:
            built = build_references(c.polarity_pattern, "LH", strand, c.GT_resids)
        builder.assert_called_once_with(c.polarity_pattern, "LH", layer_dir=None)
        layers = data_path("reference_templates", "layers")
        self.assertEqual([call.args[0] for call in reader.call_args_list],
                         [layers / "a.pdb", layers / "s.pdb"])
        template = build_reference_template(c.polarity_pattern, "LH")
        for tract_index, resids in enumerate(c.GT_resids.values()):
            for layer_index, resid in enumerate(resids):
                number = tract_index * len(c.polarity_pattern) + layer_index + 1
                expected = template.loc[template.res_seq.eq(number)]
                actual = built.full_ref.loc[built.full_ref.res_seq.eq(resid)]
                self.assertEqual(actual.atom_name.tolist(), expected.atom_name.tolist())
                np.testing.assert_allclose(actual[["x", "y", "z"]], expected[["x", "y", "z"]])
        self.assertEqual(set(built.ref_info), {"GT", "core", "tetrads"})
        self.assertEqual(built.full_ref_label, "rm8")
        self.assertTrue(built.full_ref.temp_factor.eq(1.0).all())
        self.assertTrue(built.full_ref.chain_id.eq("").all())
        self.assertTrue(built.full_ref.element.eq("").all())
        pd.testing.assert_frame_equal(strand, original)

    def test_custom_layer_directory_controls_reference_coordinates(self):
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
        refs = build_references("saa", "RH", strand, c.GT_resids,
                                layer_dir=layers)
        actual = refs.full_ref.loc[refs.full_ref.res_seq.eq(c.GT_resids["GT_1"][1])]
        np.testing.assert_allclose(actual[["x", "y", "z"]],
                                   expected.loc[expected.res_seq.eq(2), ["x", "y", "z"]])

    def test_pipeline_and_cli_generate_four_tetrad_lh_references(self):
        tracts = {
            f"GT_{index + 1}": list(range(4 * index + 1, 4 * index + 5))
            for index in range(4)
        }
        settings = config(
            self.root / "api", sequence="G" * 16,
            GT_resids=tracts, polarity_pattern="sass", helicity="LH",
        )
        pipeline = FoldingPipeline(settings)
        with redirect_stdout(StringIO()):
            pipeline.prepare()
        self.assertEqual(len(pipeline.references.full_ref), 16 * 11)
        self.assertEqual(len(pipeline.references.ref_info["tetrads"]), 4)

        values = settings.to_dict()
        del values["dt_ps"]
        input_path = self.root / "input.json"
        input_path.write_text(json.dumps(values))
        work = self.root / "cli"
        with redirect_stdout(StringIO()):
            status = main([
                "prepare_nogmx", "-i", str(input_path), "-w", str(work),
            ])
        self.assertEqual(status, 0)
        restored = FoldingPipeline.from_directory(work)
        pd.testing.assert_frame_equal(
            restored.references.full_ref, read_pdb(pipeline.paths.full_ref)
        )


if __name__ == "__main__":
    unittest.main()
