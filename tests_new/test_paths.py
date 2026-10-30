"""Packaged filename defaults and their use across all pipeline stages."""
from contextlib import redirect_stdout, redirect_stderr
from dataclasses import fields
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from g4_builder import FoldingPipeline, RunPaths, SimulationConfig
from g4_builder.cli import main
from g4_builder.paths import load_file_names
from g4_builder.resources import data_path
from g4_builder.logging import configure_logging
from tests_new.test_library import FakeRunner, config


class PathTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(lambda: configure_logging(console=False))
        self.root = Path(self.temp.name)
        self.defaults = dict(load_file_names())
        self.file = self.root / "file_names.json"

    def use_defaults(self, values):
        self.file.write_text(json.dumps(values))
        original = data_path
        def resource(*parts):
            return self.file if parts == ("file_names.json",) else original(*parts)
        return patch("g4_builder.paths.data_path", side_effect=resource)

    def test_default_names_and_removed_config_option(self):
        paths = RunPaths(self.root)
        self.assertEqual(paths.em_prefix, "em")
        self.assertEqual(paths.simulation_prefix, "smd")
        self.assertEqual(paths.em_tpr.name, "em.tpr")
        self.assertEqual(paths.mdp.name, "smd.mdp")
        self.assertEqual(paths.tpr.name, "smd.tpr")
        self.assertEqual(paths.frame(2).name, "frame_2.pdb")
        self.assertNotIn("prefix", {f.name for f in fields(SimulationConfig)})
        values = config(self.root).to_dict()
        values["prefix"] = "custom"
        with self.assertRaisesRegex(ValueError, "Unknown input settings"):
            SimulationConfig.from_dict(values, working_dir=self.root)

    def test_custom_default_names_flow_through_pipeline(self):
        names = dict(self.defaults, em_tpr="relax.min.tpr", tpr="fold.prod.tpr",
            strand="start.pdb", full_ref="complete.pdb", reference="guanine_{name}.pdb",
            plumed="restraints.dat", colvar="rmsds.dat", topology="system.top",
            posre="positions.itp", conf="processed.pdb", box="boxed.pdb",
            frame="best_{index}.pdb", manifest="metadata.json", log="commands.log")
        with self.use_defaults(names):
            work = self.root / "run"
            runner = FakeRunner(work)
            pipeline = FoldingPipeline(config(work), runner)
            with redirect_stdout(StringIO()):
                pipeline.prepare()
                pipeline.minimize()
            pipeline.prepare_simulation()
            pipeline.run()
            self.assertTrue((work / "metadata.json").is_file())
            self.assertTrue((work / "fold.prod.tpr").is_file())
            self.assertTrue((work / "fold.prod.mdp").is_file())
            self.assertTrue((work / "relax.min.pdb").is_file())
            self.assertTrue((work / "fold.prod.pdb").is_file())
            self.assertTrue((work / "guanine_ref_1.pdb").is_file())
            self.assertIn("REFERENCE=guanine_ref_1.pdb", pipeline.protocol.text)
            self.assertIn("FILE=rmsds.dat", pipeline.protocol.text)
            self.assertEqual(pipeline.result.outputs.trajectories[0].name, "fold.prod.part0002.xtc")
            md_runs = [args for args, _, _ in runner.calls if args[0] == "mdrun"]
            self.assertEqual(md_runs[0][md_runs[0].index("-deffnm") + 1], "relax.min")
            self.assertEqual(md_runs[1][md_runs[1].index("-deffnm") + 1], "fold.prod")
            end = pipeline.result.protocol_steps * pipeline.config.dt_ps
            pipeline.paths.colvar.write_text(
                f"#! FIELDS time {pipeline.result.full_ref_label}\n"
                f"{end + 10} 0.1\n{end + 20} 0.2\n")
            with redirect_stdout(StringIO()):
                pipeline.extract()
            self.assertEqual(pipeline.result.frame_paths[1].name, "best_1.pdb")
            restored = FoldingPipeline.from_directory(work)
            self.assertEqual(restored.result.outputs, pipeline.result.outputs)
            self.assertEqual(restored.paths.manifest.name, "metadata.json")
            self.assertEqual(restored.references.entries[0].filename, "guanine_ref_1.pdb")

    def test_cli_finds_custom_manifest(self):
        names = dict(self.defaults, manifest="metadata.json", tpr="fold.tpr")
        example = Path(__file__).resolve().parents[1] / "examples" / "input.json"
        with self.use_defaults(names), redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            arguments = ["prepare_nogmx", "-i", str(example), "-w", str(self.root / "cli")]
            self.assertEqual(main(arguments), 0)
            self.assertEqual(main(arguments), 0)
        self.assertTrue((self.root / "cli" / "metadata.json").is_file())

    def test_invalid_default_file_definitions(self):
        cases = [dict(self.defaults, em_tpr="input.txt"),
                 dict(self.defaults, tpr="../outside.tpr"),
                 dict(self.defaults, tpr="em.tpr"),
                 dict(self.defaults, frame="single.pdb"),
                 dict(self.defaults, tpr="{unknown}.tpr")]
        missing = dict(self.defaults)
        missing.pop("manifest")
        cases.append(missing)
        for values in cases:
            with self.subTest(values=values), self.use_defaults(values), self.assertRaises(ValueError):
                RunPaths(self.root)


if __name__ == "__main__":
    unittest.main()
