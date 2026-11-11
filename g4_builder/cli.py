"""Terminal client. Run metadata is authoritative for subsequent stages."""

import argparse
from dataclasses import replace
from pathlib import Path
import subprocess
import sys

from .client_input import load_input, validate_cli_config
from .logging import configure_logging
from .paths import RunPaths, require_files
from .pipeline import FoldingPipeline


def build_parser():
    """Define the four client commands and their flags."""
    parser = argparse.ArgumentParser(
        prog="g4-builder",
        description="Prepare, run and analyze a DNA folding simulation.",
    )
    modes = parser.add_subparsers(dest="mode", required=True)
    descriptions = {
        "prepare_nogmx": (
            "Build strand, references and PLUMED protocol without GROMACS."
        ),
        "prepare": (
            "Prepare all inputs through minimization and the simulation TPR."
        ),
        "run_sim": "Run the simulation using saved run.json metadata.",
        "extract": (
            "Extract low-RMSD frames using saved run.json metadata. "
            "Bare input filenames use --working-dir; relative paths with "
            "directories use the current directory."
        ),
    }
    for mode, description in descriptions.items():
        command = modes.add_parser(
            mode, help=description, description=description
        )
        command.add_argument(
            "-w", "--working-dir", required=True, type=Path,
            help="Run directory (separate from JSON input).",
        )
        if mode.startswith("prepare"):
            command.add_argument(
                "-i", "--input", type=Path, required=True,
                help="Simulation settings JSON (no execution blocks).",
            )
            command.add_argument(
                "--overwrite", action="store_true",
                help=(
                    "Explicitly regenerate inputs; does not delete old MD "
                    "outputs."
                ),
            )
        elif mode == "extract":
            command.add_argument(
                "-trj", "--trajectories", nargs="+", type=Path,
                help=(
                    "One XTC or ordered trajectory parts. "
                    "Default: recorded run outputs."
                ),
            )
            command.add_argument(
                "-s", "--structure", type=Path,
                help=(
                    "Structure/TPR input. "
                    "Default: recorded structure or prepared TPR."
                ),
            )
            command.add_argument(
                "-colvar", "--colvar", type=Path,
                help=(
                    "COLVAR input. "
                    "Default: recorded COLVAR or working_dir/COLVAR."
                ),
            )
    return parser


def _verify_input(saved, requested):
    """Reject settings that would change the prepared model or protocol."""
    old = saved.to_dict()
    new = requested.to_dict()
    # These settings do not alter the prepared molecular model or protocol.
    allowed = {"gromacs_executable", "n_frames"}
    conflicts = [
        key for key in old if key not in allowed and old[key] != new[key]
    ]
    if conflicts:
        raise ValueError(
            f"Input conflicts with run.json for: {', '.join(conflicts)}. "
            "Use the original settings, a new directory, or --overwrite "
            "with a preparation mode."
        )


def execute(args):
    """Load or restore a run, execute the requested stage and save metadata."""
    work = args.working_dir.expanduser().resolve()
    manifest = RunPaths(work).manifest
    preparing = args.mode in ("prepare_nogmx", "prepare")
    if preparing:
        requested = load_input(args.input, work)
    else:
        requested = None

    # Only preparation accepts new settings; later stages restore run.json.
    if preparing and (not manifest.exists() or args.overwrite):
        pipeline = FoldingPipeline(requested)
        configure_logging(pipeline.paths.log)
        pipeline.prepare(overwrite=args.overwrite)
    else:
        if not manifest.is_file():
            raise FileNotFoundError(
                f"Preparation metadata not found: {manifest}. "
                "Run prepare or prepare_nogmx first."
            )
        pipeline = FoldingPipeline.from_directory(work)
        validate_cli_config(pipeline.config)
        if requested is not None:
            _verify_input(pipeline.config, requested)
            pipeline.config = replace(
                pipeline.config,
                gromacs_executable=requested.gromacs_executable,
                n_frames=requested.n_frames,
            )
            pipeline.runner.executable = requested.gromacs_executable
        configure_logging(pipeline.paths.log)
    if args.mode == "prepare":
        if "minimized" not in pipeline.result.completed_stages:
            pipeline.minimize()
        else:
            require_files(pipeline.paths.em_structure, pipeline.paths.topology)
        if "simulation_prepared" not in pipeline.result.completed_stages:
            pipeline.prepare_simulation()
        else:
            require_files(pipeline.paths.tpr)
    elif args.mode == "run_sim":
        pipeline.run()
    elif args.mode == "extract":
        pipeline.extract(
            trajectories=args.trajectories,
            structure=args.structure,
            colvar=args.colvar,
        )
    # Also persist configuration when an already completed preparation was
    # reused, without rerunning any external command.
    if preparing:
        pipeline.save_metadata()
    print(f"{args.mode} completed. Run metadata: {pipeline.paths.manifest}")


def main(argv=None):
    """Run the client and return an exit code for its outcome."""
    args = build_parser().parse_args(argv)
    try:
        execute(args)
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except (
        OSError, ValueError, TypeError, KeyError, subprocess.SubprocessError
    ) as exc:
        print(f"g4-builder: {exc}", file=sys.stderr)
        return 1
    return 0
