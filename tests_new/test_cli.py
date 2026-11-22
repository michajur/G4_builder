"""CLI JSON boundaries, extraction flags and staged invocation tests."""
from contextlib import redirect_stdout, redirect_stderr
from dataclasses import fields, asdict, replace
from io import StringIO
import json
from os.path import relpath
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from g4_builder.cli import main
from g4_builder.client_input import CLI_FIXED_CONFIG, load_input
from g4_builder.config import ScheduleConfig, SimulationConfig
from g4_builder.logging import configure_logging
from tests_new.test_library import FakeRunner


class MenuRunner(FakeRunner):
    def probe(self, arguments, timeout=2):
        return "Select the Force Field:\n 4: Amber BSC1\n"


class ClientTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(lambda: configure_logging(console=False))
        self.root = Path(self.temp.name)
        self.work = self.root / "run"
        self.input = self.root / "input.json"
        self.values = {
            "sequence": "GGGTAGGGCGGCGGGGCAGGGT",
            "GT_resids": {"GT_1": [11, 1, 2], "GT_2": [8, 7, 6],
                          "GT_3": [19, 20, 21], "GT_4": [15, 14, 13]},
            "polarity_pattern": "ssa", "helicity": "RH", "force_field": 4,
            "n_frames": 2,
        }
        self.write_input()

    def write_input(self):
        self.input.write_text(json.dumps(self.values))

    def invoke(self, mode, *extra):
        arguments = [mode, "-w", str(self.work)]
        if mode.startswith("prepare"):
            arguments += ["-i", str(self.input)]
        arguments += list(extra)
        out, err = StringIO(), StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            status = main(arguments)
        return status, out.getvalue(), err.getvalue()

    def test_defaults_and_complete_example(self):
        configuration = load_input(self.input, self.work)
        self.assertEqual(configuration.dielectric_const, 15)
        self.assertEqual(configuration.saving_stride, 5000)
        self.assertEqual(configuration.after_folding_steps, 1500000)
        self.assertEqual(configuration.gromacs_executable, "gmx")
        self.assertEqual(configuration.force_field, 4)
        self.values["schedule"] = {"core_steps": 4500000}
        self.write_input()
        configuration = load_input(self.input, self.work)
        self.assertEqual(configuration.schedule.core_spacing, 3000000)
        self.assertEqual(configuration.schedule.core_steps, 4500000)
        self.assertEqual(configuration.schedule.core_start, 1000000)
        example = Path(__file__).resolve().parents[1] / "examples" / "input.json"
        configuration = load_input(example, self.work)
        raw = json.loads(example.read_text())
        self.assertEqual(set(raw), {f.name for f in fields(SimulationConfig)} -
                         {"working_dir"} - set(CLI_FIXED_CONFIG))
        self.assertEqual(raw["schedule"].keys(), asdict(ScheduleConfig()).keys())
        self.assertEqual(configuration.saving_stride, 5000)

    def test_only_simulation_fields_are_accepted(self):
        for key in ("working_dir", "unexpected", "run", "preparation", "extract", "logging"):
            self.input.write_text(json.dumps(dict(self.values, **{key: {}})))
            with self.subTest(key=key), self.assertRaises(ValueError):
                load_input(self.input, self.work)
        self.input.write_text('{"force_field": 1, "force_field": 2}')
        with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
            load_input(self.input, self.work)

    def test_offline_preparation_and_module_entry_point(self):
        with patch("g4_builder.runner.subprocess.Popen") as popen:
            status, output, error = self.invoke("prepare_nogmx")
        self.assertEqual(status, 0, error)
        popen.assert_not_called()
        self.assertTrue((self.work / "run.json").is_file())
        self.assertFalse((self.work / "smd.tpr").exists())
        self.assertIn("prepare_nogmx completed", output)
        result = subprocess.run([sys.executable, "-m", "g4_builder", "run_sim",
                                 "-w", str(self.work)], text=True, capture_output=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("prepare_simulation", result.stderr)

    def test_modes_reuse_simulation_metadata(self):
        runners = []
        def create_runner(directory, executable):
            runner = MenuRunner(directory)
            runners.append(runner)
            return runner
        with patch("g4_builder.pipeline.GromacsRunner", side_effect=create_runner):
            self.assertEqual(self.invoke("prepare_nogmx")[0], 0)
            status, _, error = self.invoke("prepare")
            self.assertEqual(status, 0, error)
            self.assertTrue((self.work / "smd.tpr").is_file())
            self.assertEqual(runners[-1].calls[0][2], 60)
            self.assertEqual(self.invoke("prepare")[0], 0)
            self.assertEqual(runners[-1].calls, [])
            status, _, error = self.invoke("run_sim")
            self.assertEqual(status, 0, error)
            self.assertIn("-cpi", runners[-1].calls[0][0])
            metadata = json.loads((self.work / "run.json").read_text())
            end = metadata["protocol_steps"] * metadata["config"]["dt_ps"]
            (self.work / "COLVAR").write_text(
                f"#! FIELDS time {metadata['full_ref_label']}\n"
                f"{end + 10} 0.3\n{end + 20} 0.1\n")
            status, _, error = self.invoke("extract")
            self.assertEqual(status, 0, error)
            self.assertTrue((self.work / "frame_1.pdb").exists())
            self.assertEqual(runners[-1].calls[0][1], "System\n")
            self.assertIsNone(runners[-1].calls[0][2])
            metadata = json.loads((self.work / "run.json").read_text())
            self.assertTrue((self.work / "gromacs.log").is_file())
            self.assertEqual(metadata["completed_stages"],
                             ["prepared", "minimized", "simulation_prepared", "ran", "extracted"])

    def test_extract_file_flags_single_and_multiple_trajectories(self):
        self.assertEqual(self.invoke("prepare_nogmx")[0], 0)
        directory = self.work / "trjs"
        directory.mkdir()
        for name in ("part1.xtc", "part2.xtc", "structure.pdb"):
            (directory / name).touch()
        metadata = json.loads((self.work / "run.json").read_text())
        end = metadata["protocol_steps"] * metadata["config"]["dt_ps"]
        (directory / "rmsd.dat").write_text(
            f"#! FIELDS time {metadata['full_ref_label']}\n"
            f"{end + 10} 0.3\n{end + 20} 0.1\n")
        runners = []
        def create_runner(directory, executable):
            runner = MenuRunner(directory)
            runners.append(runner)
            return runner
        with patch("g4_builder.pipeline.GromacsRunner", side_effect=create_runner):
            for input_dir in (directory, Path(relpath(directory))):
                for names in (["part1.xtc"], ["part1.xtc", "part2.xtc"]):
                    trajectories = [str(input_dir / name) for name in names]
                    with self.subTest(trajectories=trajectories):
                        status, _, error = self.invoke("extract", "-trj", *trajectories,
                                                      "-s", str(input_dir / "structure.pdb"),
                                                      "-colvar", str(input_dir / "rmsd.dat"))
                        self.assertEqual(status, 0, error)
                        calls = runners[-1].calls
                        self.assertEqual(calls[0][0][0], "trjcat" if len(trajectories) > 1 else "trjconv")
                        extraction = next(args for args, _, _ in calls if args[0] == "trjconv")
                        self.assertEqual(extraction[extraction.index("-s") + 1], str(directory / "structure.pdb"))
                        first_time = str(int(end + 20))
                        self.assertEqual(extraction[extraction.index("-dump") + 1], first_time)

    def test_input_json_is_only_a_preparation_argument(self):
        for mode in ("run_sim", "extract"):
            with self.subTest(mode=mode), redirect_stderr(StringIO()), self.assertRaises(SystemExit) as caught:
                self.invoke(mode, "-i", str(self.input))
            self.assertEqual(caught.exception.code, 2)
        with redirect_stderr(StringIO()), self.assertRaises(SystemExit) as caught:
            self.invoke("extract", "-trj")
        self.assertEqual(caught.exception.code, 2)

    def test_missing_metadata_and_preparation_conflicts(self):
        status, _, error = self.invoke("extract")
        self.assertEqual(status, 1)
        self.assertIn("Preparation metadata not found", error)
        self.assertEqual(self.invoke("prepare_nogmx")[0], 0)
        status, _, error = self.invoke("extract")
        self.assertEqual(status, 1)
        self.assertIn("supply trajectories", error)
        self.values["saving_stride"] = 200
        self.write_input()
        status, _, error = self.invoke("prepare")
        self.assertEqual(status, 1)
        self.assertIn("conflicts with run.json", error)

    def test_preparation_can_update_executable_and_frame_count(self):
        self.assertEqual(self.invoke("prepare_nogmx")[0], 0)
        self.values.update(gromacs_executable="alternate-gmx", n_frames=3)
        self.write_input()
        status, _, error = self.invoke("prepare_nogmx")
        self.assertEqual(status, 0, error)
        metadata = json.loads((self.work / "run.json").read_text())
        self.assertEqual(metadata["config"]["gromacs_executable"], "alternate-gmx")
        self.assertEqual(metadata["config"]["n_frames"], 3)

    def test_cli_fixed_fields_and_api_configurability(self):
        for name, value in CLI_FIXED_CONFIG.items():
            self.input.write_text(json.dumps(dict(self.values, **{name: value})))
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "fixed for CLI"):
                load_input(self.input, self.work)
        self.write_input()
        configuration = load_input(self.input, self.work)
        changed = replace(configuration, dt_ps=0.004)
        self.assertEqual(changed.dt_ps, 0.004)

    def test_cli_rejects_nondefault_api_model_in_saved_run(self):
        self.assertEqual(self.invoke("prepare_nogmx")[0], 0)
        manifest = self.work / "run.json"
        metadata = json.loads(manifest.read_text())
        metadata["config"]["dt_ps"] = 0.004
        manifest.write_text(json.dumps(metadata))
        status, _, error = self.invoke("run_sim")
        self.assertEqual(status, 1)
        self.assertIn("API-only settings", error)


if __name__ == "__main__":
    unittest.main()
