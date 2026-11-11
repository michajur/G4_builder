"""Select and extract low-RMSD post-protocol frames."""

from pathlib import Path

import numpy as np
import pandas as pd

from .config import integer, positive
from .paths import require_files


class TrajectoryFrameExtractor:
    """Select frames, then extract them using the shared runner and its paths.

    Only the runner is configured at construction. Pass COLVAR selection inputs
    to find_best_frames(), then trajectory and structure inputs to
    extract_frames().
    The selected times/RMSDs and saved paths remain available on the object.

    COLVAR times are ps. Multiple trajectory parts are joined with trjcat
    before extraction. COLVAR must cover those trajectories. trjconv selects
    the saved frame nearest each requested integer time; COLVAR stride should
    match the trajectory stride to make their samples correspond.

    Bare input filenames resolve from the runner's working directory. Relative
    paths containing a directory resolve from the caller's current directory.
    Absolute paths and user-home expansion are supported for every input.
    """

    def __init__(self, runner):
        self.runner = runner
        self.paths = runner.paths
        self.best_frames = {}
        self.saved_paths = {}

    @property
    def working_dir(self):
        return self.paths.working_dir

    def _resolve_path(self, path):
        """Interpret only bare filenames relative to the run directory."""
        path = Path(path).expanduser()
        if not path.is_absolute() and path.parent == Path("."):
            path = self.working_dir / path
        return path.resolve()

    def find_best_frames(
        self, rmsd_label, protocol_steps, colvar=None, n=5, dt_ps=0.002
    ):
        """Read COLVAR and return/store the n lowest-RMSD post-protocol times.

        No trajectories, structure file or GROMACS executable are needed here.
        colvar defaults to paths.colvar. Each call clears previous selection
        and saved-path results, including when validation or selection fails.
        """
        self.best_frames = {}
        self.saved_paths = {}
        protocol_steps = integer("protocol_steps", protocol_steps)
        dt_ps = positive("dt_ps", dt_ps)
        n = integer("n", n, 1)
        protocol_end_ps = protocol_steps * dt_ps
        if colvar is None:
            colvar = self.paths.colvar
        else:
            colvar = self._resolve_path(colvar)
        with colvar.open() as handle:
            header = handle.readline().strip()
        labels = header.removeprefix("#! FIELDS").lstrip("# ").split()
        if (
            len(labels) < 2
            or len(labels) != len(set(labels))
            or rmsd_label not in labels[1:]
        ):
            raise ValueError(
                "COLVAR header must contain unique fields including "
                f"{rmsd_label!r}"
            )
        frame = pd.read_csv(
            colvar, sep=r"\s+", skiprows=1, comment="#", names=labels
        )
        times = pd.to_numeric(
            frame.iloc[:, 0], errors="raise"
        ).to_numpy(dtype=float)
        rmsds = pd.to_numeric(
            frame[rmsd_label], errors="raise"
        ).to_numpy(dtype=float)

        # Filter before rounding, then check the boundary again: rounding must
        # not move an otherwise eligible sample back into the folding protocol.
        mask = (
            np.isfinite(times)
            & np.isfinite(rmsds)
            & (times >= protocol_end_ps)
        )
        candidates = pd.DataFrame({
            "time": [int(round(time)) for time in times[mask]],
            "rmsd": rmsds[mask],
        })
        candidates = candidates.loc[candidates.time >= protocol_end_ps]
        # Sort first to keep the lowest RMSD at each rounded time.
        candidates = candidates.sort_values("rmsd", kind="stable")
        candidates = candidates.drop_duplicates("time")
        if len(candidates) < n:
            raise ValueError(
                f"Requested {n} frames; only {len(candidates)} distinct "
                "post-protocol times remain"
            )
        self.best_frames = {
            int(row.time): float(row.rmsd)
            for row in candidates.head(n).itertuples(index=False)
        }
        return self.best_frames

    def extract_frames(self, trajectories, structure, timeout=None):
        """Extract the System group at the previously selected times.

        trajectories is one path or an ordered iterable of parts. structure is
        the matching TPR/structure path. Results are printed and stored in
        saved_paths; the method returns None.
        """
        self.saved_paths = {}
        if not self.best_frames:
            raise ValueError("Call find_best_frames() first")
        if isinstance(trajectories, (str, Path)):
            trajectories = [trajectories]
        trajectories = tuple(self._resolve_path(path) for path in trajectories)
        if not trajectories:
            raise ValueError("At least one trajectory is required")
        structure = self._resolve_path(structure)
        require_files(structure, *trajectories)
        trajectory = trajectories[0]
        if len(trajectories) > 1:
            trajectory = self.paths.joined_trajectory
            if trajectory in trajectories:
                raise ValueError(
                    "The joined trajectory output must not replace an input"
                )
            self.runner.run(
                ["trjcat", "-f", *trajectories, "-o", trajectory],
                timeout=timeout,
            )
            require_files(trajectory)
        for index, time in enumerate(self.best_frames, 1):
            path = self.paths.frame(index)
            self.runner.run(
                ["trjconv", "-f", trajectory, "-s", structure,
                 "-dump", time, "-o", path],
                input_text="System\n",
                timeout=timeout,
            )
            require_files(path)
            self.saved_paths[index] = path
            print(f"Saved frame_{index}: {path}")
