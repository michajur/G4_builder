"""Tests for the new package; legacy modules and notebooks are not imported."""
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
import json
import logging
from os.path import relpath
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from g4_builder import (SimulationConfig, ScheduleConfig, FoldingPipeline,
    DNAStrandBuilder, build_references, ReferenceSet, build_protocol, read_pdb, save_pdb,
    GromacsRunner, GromacsPreparation, TrajectoryFrameExtractor, configure_logging)


def config(directory, **changes):
    values = dict(sequence="GGGTAGGGCGGCGGGGCAGGGT",
                  GT_resids={"GT_1": [11, 1, 2], "GT_2": [8, 7, 6],
                             "GT_3": [19, 20, 21], "GT_4": [15, 14, 13]},
                  polarity_pattern="saa", helicity="RH", working_dir=directory,
                  force_field="amber99bsc1", n_frames=2)
    values.update(changes)
    return SimulationConfig(**values)


class FakeRunner(GromacsRunner):
    """Act like GROMACS at the file boundary, with changed atom serials."""
    def __init__(self, working_dir, reorder_atoms=True):
        super().__init__(working_dir)
        self.calls = []
        self.reorder_atoms = reorder_atoms

    def run(self, arguments, input_text="", timeout=None):
        args = list(map(str, arguments))
        self.calls.append((args, input_text, timeout))
        def value(flag): return args[args.index(flag) + 1]
        if args[0] == "pdb2gmx":
            strand = read_pdb(value("-f"))
            # Exercise matching by identity, not atom order or old serial.
            if self.reorder_atoms:
                strand = strand.iloc[::-1].reset_index(drop=True)
                strand["serial"] = np.arange(1, len(strand) + 1)
            save_pdb(strand, value("-o"))
            Path(value("-p")).touch()
        elif args[0] == "editconf":
            shutil.copyfile(value("-f"), value("-o"))
        elif args[0] in ("grompp", "trjconv", "trjcat"):
            Path(value("-o")).touch()
        elif args[0] == "mdrun":
            if "-plumed" not in args:
                output = self.working_dir / value("-c")
                output.touch()
            else:
                trajectory = self.working_dir / f"{value('-deffnm')}.xtc"
                trajectory.touch()
                Path(value("-c")).touch()
                self.paths.colvar.touch()


class LibraryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def quiet_prepare(self, pipeline):
        with redirect_stdout(StringIO()):
            return pipeline.prepare()

    def test_components_share_runner_paths_without_running_commands(self):
        work = self.root / "not_created"
        runner = GromacsRunner(work, executable="unavailable-gromacs")
        preparation = GromacsPreparation(runner)
        extractor = TrajectoryFrameExtractor(runner)
        pipeline = FoldingPipeline(config(work, dt_ps=0.004), runner)
        for component in (preparation, extractor, pipeline, pipeline.gromacs):
            self.assertIs(component.paths, runner.paths)
        self.assertEqual(extractor.working_dir, runner.working_dir)
        self.assertEqual(extractor.best_frames, {})
        self.assertEqual(extractor.saved_paths, {})
        self.assertIsNone(preparation.dt_ps)
        self.assertFalse(work.exists())

    def test_pipeline_rejects_runner_in_another_directory(self):
        work = self.root / "run"
        runner = GromacsRunner(self.root / "other")
        with self.assertRaisesRegex(ValueError, "working directories must match"):
            FoldingPipeline(config(work), runner)
        self.assertFalse(work.exists())
        self.assertFalse(runner.working_dir.exists())

    def test_run_records_single_trajectory_and_ignores_old_parts(self):
        runner = FakeRunner(self.root)
        preparation = GromacsPreparation(runner)
        paths = runner.paths
        paths.tpr.touch()
        paths.plumed.touch()
        prefix = paths.simulation_prefix
        for suffix in (".part0010.xtc", ".part0001.xtc"):
            (self.root / f"{prefix}{suffix}").touch()
        preparation.run_sim()
        self.assertEqual(preparation.outputs.trajectories, (paths.trajectory,))
        self.assertEqual(paths.trajectory.name, "smd.xtc")
        self.assertNotIn("-noappend", runner.calls[-1][0])
        self.assertIn("-cpi", runner.calls[-1][0])

    def test_run_requires_single_trajectory_and_colvar(self):
        runner = FakeRunner(self.root)
        preparation = GromacsPreparation(runner)
        paths = runner.paths
        paths.tpr.touch()
        paths.plumed.touch()
        # Old trajectory parts must not conceal a missing new output.
        (self.root / "smd.part0001.xtc").touch()
        for missing in (paths.trajectory, paths.colvar):
            with self.subTest(missing=missing):
                paths.trajectory.touch()
                paths.colvar.touch()
                missing.unlink()
                with patch.object(runner, "run"), self.assertRaises(FileNotFoundError):
                    preparation.run_sim()
                self.assertIsNone(preparation.outputs)

    def test_reference_set_saves_full_and_partial_references(self):
        c = config(self.root)
        strand = DNAStrandBuilder().build(c.sequence)
        references = build_references(c.polarity_pattern, c.helicity, strand, c.GT_resids)
        self.assertIsInstance(references, ReferenceSet)
        references_dir = self.root / "references"
        full_path = self.root / "full_ref.pdb"
        with redirect_stdout(StringIO()) as output:
            self.assertIsNone(references.save_full_ref(full_path))
            self.assertIsNone(references.save_partial_refs(references_dir))
        self.assertEqual(set(references.saved_paths), set(references.partial_refs))
        for entry in references.entries:
            references_path = references.saved_paths[entry.name]
            self.assertEqual(references_path, references_dir / entry.filename)
            pd.testing.assert_frame_equal(read_pdb(references_path), entry.atoms)
            self.assertIn(str(references_path), output.getvalue())
        for path in [full_path, *references.saved_paths.values()]:
            text = path.read_text()
            self.assertNotIn("TER", text)
            self.assertNotIn("END", text)
            atoms = read_pdb(path)
            self.assertTrue(atoms.temp_factor.eq(1.0).all())
            self.assertTrue(atoms.chain_id.eq("").all())
            self.assertTrue(atoms.element.eq("").all())
        pd.testing.assert_frame_equal(read_pdb(full_path), references.full_ref)

    def test_validation_precedes_output_and_copies_inputs(self):
        destination = self.root / "not_created"
        bad = [dict(sequence="GGNX"), dict(GT_resids={"a": [1, 2]}),
               dict(polarity_pattern="xx"), dict(helicity="right"),
               dict(force_field=0), dict(saving_stride=True), dict(dt_ps=0),
               dict(dielectric_const=float("nan"))]
        for changes in bad:
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                config(destination, **changes)
            self.assertFalse(destination.exists())
        tracts = {"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]}
        c = config(destination, sequence="G" * 8, polarity_pattern="ss", GT_resids=tracts)
        tracts["a"][0] = 99
        self.assertEqual(c.GT_resids["a"], (1, 2))
        pipeline = FoldingPipeline(c)
        self.assertFalse(destination.exists())
        with patch("g4_builder.runner.subprocess.Popen") as popen:
            self.quiet_prepare(pipeline)
        popen.assert_not_called()

    def test_offline_prepare_metadata_and_fresh_session(self):
        p = FoldingPipeline(config(self.root / "run"))
        result = self.quiet_prepare(p)
        self.assertEqual(result.protocol_steps, 12000000)
        self.assertEqual(result.simulation_steps, 14500000)
        self.assertEqual(result.full_ref_label, "rm8")
        self.assertEqual(len(p.references.entries), 11)
        for entry in p.references.entries:
            atoms = read_pdb(result.reference_paths[entry.name])
            self.assertTrue(atoms.temp_factor.eq(1.0).all())
            self.assertTrue(atoms.chain_id.eq("").all())
            self.assertTrue(atoms.element.eq("").all())
            self.assertNotIn("TER", result.reference_paths[entry.name].read_text())
        resumed = FoldingPipeline.from_directory(p.paths.working_dir)
        self.assertEqual(resumed.result.protocol_steps, result.protocol_steps)
        self.assertEqual(resumed.config.to_dict(), p.config.to_dict())
        self.assertEqual(resumed.protocol.text, p.protocol.text)
        self.assertEqual(resumed.references.saved_paths, result.reference_paths)
        with self.assertRaises(FileExistsError):
            p.prepare()
        with self.assertRaisesRegex(ValueError, "minimize"):
            p.prepare_simulation()
        with self.assertRaisesRegex(ValueError, "prepare_simulation"):
            p.run()

    def test_reference_identity_and_build_input_immutability(self):
        from g4_builder.resources import data_path
        # With two tetrads and no bulges, GT + core + tetrads yields nine refs.
        c = config(self.root, sequence="G" * 8, polarity_pattern="ss",
                   GT_resids={"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]})
        strand = DNAStrandBuilder().build(c.sequence)
        original_strand = strand.copy(deep=True)
        template = read_pdb(data_path("reference_templates", "RH_2tet", "ss.pdb"))
        original_template = template.copy(deep=True)
        tracts = {name: list(values) for name, values in c.GT_resids.items()}
        with patch("g4_builder.references.read_pdb", return_value=template):
            references = build_references("ss", "RH", strand, tracts)
        self.assertEqual(references.full_ref_label, "rm7")
        self.assertEqual(sum(e.is_full_core for e in references.entries), 1)
        pd.testing.assert_frame_equal(template, original_template)
        pd.testing.assert_frame_equal(strand, original_strand)
        self.assertEqual(tracts["a"], [1, 2])
        tracts["a"][0] = 99
        self.assertEqual(references.entries[0].resids, (1, 2))
        # Duplicating an equal-sized table does not alter explicit identity.
        from g4_builder.references import Reference
        first = references.entries[0]
        duplicate = Reference("extra", "GT", first.resids, "other", references.full_ref.copy())
        references.entries = (duplicate,) + references.entries
        self.assertEqual(references.full_ref_label, "rm7")

    def test_missing_template_and_atom_mismatch(self):
        from g4_builder.resources import data_path
        c = config(self.root)
        strand = DNAStrandBuilder().build(c.sequence)
        with self.assertRaises(FileNotFoundError):
            build_references("saa", "RH", strand, c.GT_resids, template_dir=self.root)
        template = read_pdb(data_path("reference_templates", "RH_3tet", "saa.pdb"))
        template.loc[0, "atom_name"] = "WRONG"
        with patch("g4_builder.references.read_pdb", return_value=template):
            with self.assertRaisesRegex(ValueError, "atom names"):
                build_references("saa", "RH", strand, c.GT_resids)

    def test_protocol_counts_offsets_and_custom_schedule(self):
        p = FoldingPipeline(config(self.root / "run"))
        self.quiet_prepare(p)
        text = p.protocol.text
        self.assertIn("RMSD REFERENCE=ref_1.pdb TYPE=OPTIMAL LABEL=rm1", text)
        self.assertEqual(text.count("RMSD REFERENCE="), len(p.references.entries))
        self.assertEqual(text.count("MOVINGRESTRAINT ..."), 6)
        self.assertIn("STEP0=1000000  AT0=0.5", text)
        self.assertIn("STEP0=3000000  AT0=0.5", text)
        self.assertIn("STEP0=9000000  AT0=0.1", text)
        # Generation is pure; the schedule is independent of writing paths.
        schedule = ScheduleConfig(core_spacing=3000000)
        changed = build_protocol(p.references, p.strand.serial, schedule)
        self.assertEqual(changed.final_steps, 16000000)
        with self.assertRaises(ValueError):
            ScheduleConfig(core=((2, 0, 1), (1, 0, 0)))

    def test_pipeline_stages_and_recorded_outputs(self):
        work = self.root / "run"
        runner = FakeRunner(work)
        p = FoldingPipeline(config(work, dt_ps=0.004, dielectric_const=80, saving_stride=100), runner)
        self.quiet_prepare(p)
        with redirect_stdout(StringIO()) as output, patch.object(
                p.references, "remap_serials", wraps=p.references.remap_serials) as remap:
            p.minimize()
        remap.assert_called_once()
        self.assertIn("changed after pdb2gmx; remapping", output.getvalue())
        self.assertFalse(p.references.needs_serial_remap(read_pdb(p.paths.conf)))
        processed = read_pdb(p.paths.conf).set_index(["res_seq", "atom_name"])
        for atom in p.references.full_ref.itertuples():
            self.assertEqual(atom.serial, processed.loc[(atom.res_seq, atom.atom_name), "serial"])
        p.prepare_simulation()
        mdp = p.paths.mdp.read_text()
        self.assertIn("epsilon_r               = 80", mdp)
        self.assertRegex(mdp, r"dt\s*=\s*0.004")
        self.assertIn("nsteps                  = 14500000", mdp)
        p.run()
        self.assertEqual(p.result.outputs.trajectories[0].name, "smd.xtc")
        self.assertEqual(p.result.outputs.structure, p.paths.tpr)
        end = p.result.protocol_steps * p.config.dt_ps
        p.paths.colvar.write_text(f"#! FIELDS time {p.result.full_ref_label}\n"
                                 f"0 0.0\n{end} 0.4\n{end + 10} 0.1\n{end + 20} 0.2\n")
        with redirect_stdout(StringIO()):
            result = p.extract()
        self.assertEqual(list(result.best_frames), [int(end + 10), int(end + 20)])
        self.assertEqual(len(result.frame_paths), 2)
        self.assertEqual(runner.calls[-1][1], "System\n")
        restored = FoldingPipeline.from_directory(work, runner)
        self.assertEqual(restored.result.frame_paths, result.frame_paths)
        self.assertEqual(restored.result.outputs, result.outputs)
        self.assertEqual(restored.result.best_frames, result.best_frames)
        destination = self.root / "moved"
        shutil.copytree(work, destination)
        moved = FoldingPipeline.from_directory(destination)
        self.assertEqual(moved.result.outputs.structure, destination / "smd.tpr")

    def test_gt_hold_kappa_follows_last_initial_stage(self):
        p = FoldingPipeline(config(self.root / "run"))
        self.quiet_prepare(p)
        for kappa in (9000, 0):
            schedule = ScheduleConfig(gt_initial=((0, 0.3, 0), (400000, 0.0, kappa)))
            self.assertEqual(schedule.gt_hold_kappa, kappa)
            protocol = build_protocol(p.references, p.strand.serial, schedule)
            gt_block = protocol.text.split("MOVINGRESTRAINT ...", 1)[1].split("... MOVINGRESTRAINT", 1)[0]
            repeated = ",".join([str(kappa)] * len(p.references.ref_info["GT"]))
            self.assertIn(f"KAPPA1={repeated}\n", gt_block)
            self.assertIn(f"KAPPA2={repeated}\n", gt_block)
            saved = replace(p.config, schedule=schedule).to_dict()
            self.assertNotIn("gt_hold_kappa", saved["schedule"])
            restored = SimulationConfig.from_dict(saved, working_dir=self.root)
            self.assertEqual(restored.schedule.gt_hold_kappa, kappa)

    def test_missing_stage_files_and_explicit_extraction(self):
        p = FoldingPipeline(config(self.root / "run"), FakeRunner(self.root / "run"))
        self.quiet_prepare(p)
        with self.assertRaisesRegex(ValueError, "supply trajectories"):
            p.extract()
        with self.assertRaises(FileNotFoundError):
            p.gromacs.prep_sim_file(1000)

    def test_api_execution_arguments_remain_configurable(self):
        work = self.root / "api_controls"
        runner = FakeRunner(work)
        p = FoldingPipeline(config(work, dt_ps=0.004), runner)
        self.quiet_prepare(p)
        with redirect_stdout(StringIO()):
            p.minimize(pdb2gmx_timeout=12, discovery_timeout=3)
        self.assertEqual(runner.calls[0][2], 12)
        p.prepare_simulation()
        em_command = next(args for args, _, _ in runner.calls if args[0] == "grompp")
        self.assertEqual(em_command[-2:], ["-maxwarn", "1"])
        self.assertEqual(runner.calls[-1][0][-2:], ["-maxwarn", "2"])
        p.run(timeout=15)
        self.assertIn("-cpi", runner.calls[-1][0])
        self.assertEqual(runner.calls[-1][2], 15)
        end = p.result.protocol_steps * p.config.dt_ps
        p.paths.colvar.write_text(f"#! FIELDS time {p.result.full_ref_label}\n{end + 10} 0.1\n")
        with redirect_stdout(StringIO()):
            p.extract(n=1, timeout=15)
        self.assertEqual(runner.calls[-1][1], "System\n")
        self.assertEqual(runner.calls[-1][2], 15)

    def test_mdp_timestep_token_and_default(self):
        from g4_builder.resources import data_path
        work = self.root / "mdp_tokens"
        runner = FakeRunner(work)
        p = FoldingPipeline(config(work), runner)
        self.quiet_prepare(p)
        p.paths.topology.touch()
        p.paths.em_structure.touch()
        template = data_path("mdp", "SMD.mdp")
        original = template.read_text()
        self.assertIn("DT_PS", original)
        p.gromacs.prep_sim_file(1000)
        rendered = p.paths.mdp.read_text()
        self.assertRegex(rendered, r"dt\s*=\s*0.002")
        self.assertNotIn("DT_PS", rendered)
        self.assertEqual(template.read_text(), original)
        invalid = self.root / "invalid.mdp"
        invalid.write_text(original.replace("DT_PS", "0.002"))
        runner.calls.clear()
        with patch("g4_builder.paths.data_path", return_value=invalid):
            with self.assertRaisesRegex(ValueError, "Missing MDP template token: DT_PS"):
                p.gromacs.prep_sim_file(1000)
        self.assertEqual(runner.calls, [])

    def test_gt_release_matches_second_tetrad_stage(self):
        p = FoldingPipeline(config(self.root / "release"))
        self.quiet_prepare(p)
        schedule = ScheduleConfig(tetrads=((0, 0.1, 0), (1700000, 0.0, 20000),
                                           (4000000, 0.0, 5000)))
        self.assertEqual(schedule.gt_release_offset, 1700000)
        text = build_protocol(p.references, p.strand.serial, schedule).text
        blocks = text.split("MOVINGRESTRAINT ...")[1:]
        start = schedule.core_start + len(p.references.ref_info["core"]) * schedule.core_spacing
        held = ",".join([str(schedule.gt_hold_kappa)] * len(p.references.ref_info["GT"]))
        self.assertIn(f"STEP3={start}  AT3=", blocks[0])
        self.assertIn(f"KAPPA3={held}\n", blocks[0])
        self.assertIn(f"STEP4={start + 1700000}  AT4=", blocks[0])
        self.assertIn("KAPPA4=" + ",".join(["0"] * len(p.references.ref_info["GT"])) + "\n", blocks[0])
        self.assertIn(f"STEP1={start + 1700000}  AT1=", blocks[-1])
        saved = replace(p.config, schedule=schedule).to_dict()
        self.assertNotIn("gt_release_offset", saved["schedule"])
        self.assertEqual(SimulationConfig.from_dict(saved, working_dir=self.root).schedule.gt_release_offset, 1700000)
        with self.assertRaisesRegex(ValueError, "at least two stages"):
            ScheduleConfig(tetrads=((0, 0.1, 0),))

    def test_unchanged_reference_serials_skip_remapping(self):
        work = self.root / "unchanged"
        p = FoldingPipeline(config(work), FakeRunner(work, reorder_atoms=False))
        self.quiet_prepare(p)
        before = p.references.full_ref.copy(deep=True)
        with redirect_stdout(StringIO()) as output, patch.object(p.references, "remap_serials") as remap:
            p.minimize()
        remap.assert_not_called()
        self.assertNotIn("remapping", output.getvalue())
        pd.testing.assert_frame_equal(before, p.references.full_ref)
        # DataFrame row order alone does not change atom serial identities.
        processed = read_pdb(p.paths.conf)
        self.assertFalse(p.references.needs_serial_remap(processed.iloc[::-1]))
        # Validate mismatches instead of accepting them as an unchanged order.
        with self.assertRaisesRegex(ValueError, "missing a reference"):
            p.references.needs_serial_remap(processed.loc[~processed.atom_name.eq("N9")])
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            p.references.needs_serial_remap(pd.concat([processed, processed.iloc[:1]]))

    def test_strand_alignment_and_pdb_round_trip(self):
        builder = DNAStrandBuilder()
        original = builder.templates["T"].copy(deep=True)
        strand = builder.build("ATC")
        self.assertEqual(len(strand), 93)
        self.assertEqual(strand.groupby("res_seq").res_name.first().tolist(), ["DA", "DT", "DC"])
        np.testing.assert_allclose(builder._anchors(strand, 2),
                                   builder._anchors(builder.templates["A5"], 2), atol=1e-10)
        pd.testing.assert_frame_equal(original, builder.templates["T"])
        path = self.root / "ATC.pdb"
        builder.save_pdb(path)
        round_trip = read_pdb(path)
        np.testing.assert_allclose(round_trip[["x", "y", "z"]], strand[["x", "y", "z"]], atol=0.00051)

    def test_force_field_discovery_is_explicit(self):
        work = self.root / "run"
        runner = FakeRunner(work)
        p = FoldingPipeline(config(work), runner)
        self.assertEqual(runner.calls, [])
        self.quiet_prepare(p)
        runner.probe = lambda arguments, timeout: "Select the Force Field:\n 1: Amber BSC1\n 2: Other\n"
        self.assertEqual(p.gromacs.discover_force_fields(), {1: "Amber BSC1", 2: "Other"})
        p.gromacs.pdb2gmx(1)
        self.assertEqual(runner.calls[-1][1], "1\n")
        self.assertIn("tip3p", runner.calls[-1][0])

    def test_pdb2gmx_custom_strand_is_used_for_execution_and_discovery(self):
        custom = self.root / "inputs" / "custom.pdb"
        custom.parent.mkdir()
        save_pdb(DNAStrandBuilder().build("ATC"), custom)
        for force_field, strand in (("amber99bsc1", "inputs/custom.pdb"), (4, custom)):
            with self.subTest(force_field=force_field):
                runner = FakeRunner(self.root)
                preparation = GromacsPreparation(runner)
                self.assertFalse(runner.paths.strand.exists())
                with patch.object(runner, "probe", return_value="Select the Force Field:\n 4: Amber BSC1\n") as probe:
                    preparation.pdb2gmx(force_field, strand=strand)
                arguments = runner.calls[-1][0]
                self.assertEqual(arguments[arguments.index("-f") + 1], str(custom))
                self.assertTrue(runner.paths.conf.is_file())
                self.assertTrue(runner.paths.topology.is_file())
                if force_field == 4:
                    probe.assert_called_once()
                    arguments = probe.call_args.args[0]
                    self.assertEqual(arguments[arguments.index("-f") + 1], custom)
                    self.assertEqual(runner.calls[-1][1], "4\n")
                else:
                    probe.assert_not_called()

    def test_menu_number_configuration_and_automatic_discovery(self):
        work = self.root / "numbered"
        runner = FakeRunner(work)
        runner.probe = lambda arguments, timeout: "Select the Force Field:\n 4: Amber BSC1\n"
        c = config(work, force_field=np.int64(4))
        self.assertEqual(c.force_field, 4)
        self.assertIsInstance(c.force_field, int)
        p = FoldingPipeline(c, runner)
        self.quiet_prepare(p)
        with patch.object(runner, "probe", wraps=runner.probe) as probe:
            with redirect_stdout(StringIO()):
                p.minimize()
            probe.assert_called_once()
        self.assertEqual(runner.calls[0][1], "4\n")
        self.assertNotIn("-ff", runner.calls[0][0])
        self.assertEqual(FoldingPipeline.from_directory(work).config.force_field, 4)
        before = len(runner.calls)
        with self.assertRaisesRegex(ValueError, "not available"):
            p.gromacs.pdb2gmx(99)
        self.assertEqual(len(runner.calls), before)
        for value in (True, 0, -1, 4.5):
            with self.subTest(value=value), self.assertRaises(ValueError):
                config(work, force_field=value)

    def test_invalid_config_does_not_touch_workdir(self):
        work = self.root / "untouched"
        with self.assertRaisesRegex(ValueError, "ordered"):
            config(work, GT_resids={"a": {11, 1, 2}, "b": [8, 7, 6],
                                   "c": [19, 20, 21], "d": [15, 14, 13]})
        self.assertFalse(work.exists())

    def test_simulation_config_from_json(self):
        example = Path(__file__).resolve().parents[1] / "examples" / "input.json"
        loaded = SimulationConfig.from_json(example, working_dir=self.root / "loaded")
        self.assertEqual(loaded.polarity_pattern, "ssa")
        self.assertEqual(loaded.GT_resids["GT_1"], (11, 1, 2))
        self.assertEqual(loaded.force_field, 4)
        self.assertEqual(loaded.dielectric_const, 15)
        self.assertFalse(loaded.working_dir.exists())
        values = loaded.to_dict()
        self.assertNotIn("working_dir", values)
        values.update(dt_ps=0.004)
        path = self.root / "api.json"
        path.write_text(json.dumps(values))
        restored = SimulationConfig.from_json(path, working_dir=loaded.working_dir)
        self.assertEqual(restored.dt_ps, 0.004)
        self.assertEqual(restored.working_dir, loaded.working_dir)
        relocated = SimulationConfig.from_json(path, working_dir=self.root / "other_run")
        self.assertEqual(relocated.working_dir, self.root / "other_run")
        with self.assertRaisesRegex(TypeError, "working_dir"):
            SimulationConfig.from_json(example)
        with self.assertRaisesRegex(ValueError, "working_dir"):
            SimulationConfig.from_json(path, working_dir=None)
        path.write_text(json.dumps(dict(values, working_dir=str(self.root))))
        with self.assertRaisesRegex(ValueError, "working_dir must be supplied separately"):
            SimulationConfig.from_json(path, working_dir=self.root / "other_run")
        values["unexpected"] = 1
        path.write_text(json.dumps(values))
        with self.assertRaisesRegex(ValueError, "Unknown input settings"):
            SimulationConfig.from_json(path, working_dir=self.root)
        path.write_text('{"sequence": "GG", "sequence": "GGG"}')
        with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
            SimulationConfig.from_json(path, working_dir=self.root)

    def test_simulation_config_from_dict_requires_separate_working_directory(self):
        values = config(self.root).to_dict()
        before = json.dumps(values, sort_keys=True)
        destination = self.root / "new_run"
        loaded = SimulationConfig.from_dict(values, working_dir=destination)
        self.assertEqual(loaded.working_dir, destination)
        self.assertEqual(json.dumps(values, sort_keys=True), before)
        self.assertEqual(loaded.to_dict(), values)
        self.assertFalse(destination.exists())
        with self.assertRaisesRegex(TypeError, "working_dir"):
            SimulationConfig.from_dict(values)
        with self.assertRaisesRegex(ValueError, "working_dir"):
            SimulationConfig.from_dict(values, working_dir=None)
        for embedded in (str(destination), str(self.root / "different_run")):
            with self.subTest(embedded=embedded), self.assertRaisesRegex(
                    ValueError, "working_dir must be supplied separately"):
                SimulationConfig.from_dict(dict(values, working_dir=embedded), working_dir=destination)

    def test_config_accepts_only_current_fields(self):
        configuration = config(self.root, dielectric_const=80, saving_stride=100)
        saved = configuration.to_dict()
        self.assertIn("dielectric_const", saved)
        self.assertIn("saving_stride", saved)
        restored = SimulationConfig.from_dict(saved, working_dir=self.root)
        self.assertEqual(restored.to_dict(), configuration.to_dict())
        for name in ("unexpected", "preparation", "run", "extract", "logging"):
            invalid = dict(saved, **{name: {}})
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Unknown input settings"):
                SimulationConfig.from_dict(invalid, working_dir=self.root)
            path = self.root / "invalid.json"
            path.write_text(json.dumps(invalid))
            with self.assertRaisesRegex(ValueError, "Unknown input settings"):
                SimulationConfig.from_json(path, working_dir=self.root)
        saved["schedule"]["unexpected"] = 100
        with self.assertRaisesRegex(ValueError, "Invalid input configuration"):
            SimulationConfig.from_dict(saved, working_dir=self.root)

    def test_default_schedule_after_module_reload(self):
        # Run in isolation so reloading config does not alter other tests'
        # imported dataclass identities. Retain the package's old export,
        # reproducing a notebook with stale imports after editing the module.
        values = config(self.root).to_dict()
        values.pop("schedule")
        code = (
            "import importlib\n"
            "import g4_builder as g4\n"
            "module = importlib.import_module('g4_builder.config')\n"
            "importlib.reload(module)\n"
            f"configuration = g4.SimulationConfig(working_dir={str(self.root)!r}, **{values!r})\n"
            "assert isinstance(configuration.schedule, module.ScheduleConfig)\n"
            "assert configuration.schedule.core_start == 1000000\n"
        )
        result = subprocess.run([sys.executable, "-c", code], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_extraction_accepts_absolute_relative_and_filename_inputs(self):
        work = self.root / "testing" / "manual_run"
        work.mkdir(parents=True)
        (work / "trajectory.xtc").touch()
        (work / "structure.tpr").touch()
        (work / "COLVAR").write_text("#! FIELDS time rm1\n20 0.1\n")
        cases = [
            tuple(str(work / name) for name in ("trajectory.xtc", "structure.tpr", "COLVAR")),
            tuple(Path(relpath(work / name)) for name in ("trajectory.xtc", "structure.tpr", "COLVAR")),
            ("trajectory.xtc", "structure.tpr", "COLVAR"),
        ]
        for trajectory, structure, colvar in cases:
            with self.subTest(trajectory=trajectory):
                runner = FakeRunner(work)
                extractor = TrajectoryFrameExtractor(runner)
                extractor.find_best_frames("rm1", 0, colvar=colvar, n=1)
                with redirect_stdout(StringIO()):
                    extractor.extract_frames(trajectory, structure)
                self.assertEqual(extractor.saved_paths[1].parent, work)
                arguments = runner.calls[-1][0]
                self.assertEqual(arguments[arguments.index("-f") + 1], str(work / "trajectory.xtc"))
                self.assertEqual(arguments[arguments.index("-s") + 1], str(work / "structure.tpr"))

    def test_extraction_accepts_mixed_paths_in_trajectory_list(self):
        work = self.root / "run"
        work.mkdir()
        parts = [work / f"part{i}.xtc" for i in range(1, 4)]
        for path in [*parts, work / "structure.tpr"]:
            path.touch()
        (work / "COLVAR").write_text("#! FIELDS time rm1\n20 0.1\n")
        runner = FakeRunner(work)
        extractor = TrajectoryFrameExtractor(runner)
        extractor.find_best_frames("rm1", 0, n=1)
        with redirect_stdout(StringIO()):
            extractor.extract_frames(["part1.xtc", Path(relpath(parts[1])), parts[2]], "structure.tpr")
        self.assertEqual(runner.calls[0][0][0], "trjcat")
        self.assertEqual(runner.calls[0][0][2:5], list(map(str, parts)))

    def test_frame_selection_rounding_and_multiple_parts(self):
        runner = FakeRunner(self.root)
        for name in ("a.xtc", "b.xtc", "input.tpr"):
            (self.root / name).touch()
        (self.root / "COLVAR").write_text("time rm1\n9.9 0\n10 0.9\n11.2 0.2\n11.4 0.1\n12.7 0.3\n14 0.2\n")
        e = TrajectoryFrameExtractor(runner)
        self.assertEqual(list(e.find_best_frames("rm1", 2500, n=3, dt_ps=0.004).items()),
                         [(11, 0.1), (14, 0.2), (13, 0.3)])
        with redirect_stdout(StringIO()) as output:
            self.assertIsNone(e.extract_frames(["a.xtc", "b.xtc"], "input.tpr"))
        self.assertEqual(runner.calls[0][0][0], "trjcat")
        self.assertEqual(len(e.saved_paths), 3)
        self.assertIn(str(e.saved_paths[1]), output.getvalue())
        with self.assertRaisesRegex(ValueError, "distinct post-protocol"):
            e.find_best_frames("rm1", 999999)

    def test_reselection_replaces_results_and_failure_clears_old_frames(self):
        runner = FakeRunner(self.root)
        extractor = TrajectoryFrameExtractor(runner)
        (self.root / "COLVAR").write_text(
            "#! FIELDS time rm1 rm2\n10 0.2 0.3\n20 0.1 0.2\n30 0.4 0.1\n")
        # Selection works before any trajectory or structure is available.
        self.assertEqual(extractor.find_best_frames("rm1", 0, n=1), {20: 0.1})
        self.assertEqual(runner.calls, [])
        for name in ("trajectory.xtc", "input.tpr"):
            (self.root / name).touch()
        with redirect_stdout(StringIO()):
            extractor.extract_frames("trajectory.xtc", "input.tpr")
        self.assertTrue(extractor.saved_paths)
        self.assertEqual(extractor.find_best_frames("rm2", 0, n=2), {30: 0.1, 20: 0.2})
        self.assertEqual(extractor.saved_paths, {})
        with self.assertRaisesRegex(ValueError, "header"):
            extractor.find_best_frames("missing_label", 0)
        self.assertEqual(extractor.best_frames, {})
        self.assertEqual(extractor.saved_paths, {})
        calls_before = len(runner.calls)
        with self.assertRaisesRegex(ValueError, "find_best_frames"):
            extractor.extract_frames("trajectory.xtc", "input.tpr")
        self.assertEqual(len(runner.calls), calls_before)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_streams_failure_and_timeout(self):
        log_file = self.root / "commands.log"
        logger = configure_logging(log_file, console=False)
        runner = GromacsRunner(self.root, sys.executable, logger)
        self.assertIsNone(runner.run(["-c", "import sys; print('OUT'); print('ERR', file=sys.stderr)"]))
        self.assertIn("OUT", log_file.read_text())
        self.assertIn("ERR", log_file.read_text())
        with self.assertRaises(subprocess.CalledProcessError) as caught:
            runner.run(["-c", "import sys; print('BAD', file=sys.stderr); sys.exit(7)"])
        self.assertEqual(caught.exception.returncode, 7)
        self.assertIn("BAD", log_file.read_text())
        with self.assertRaises(subprocess.TimeoutExpired):
            runner.run(["-c", "import time; time.sleep(5)"], timeout=0.05)
        runner.executable = "g4-nonexistent-executable-test"
        with self.assertRaises(FileNotFoundError):
            runner.run(["mdrun"])
        for handler in list(logger.handlers):
            if getattr(handler, "_g4_owned", False):
                logger.removeHandler(handler)
                handler.close()


if __name__ == "__main__":
    unittest.main()
