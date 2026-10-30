"""Shared subprocess execution for preparation and frame extraction."""
from collections import deque
import logging
from pathlib import Path
import shutil
import subprocess
import threading

from .paths import RunPaths


class GromacsRunner:
    """Stream both GROMACS output channels to logging without returning them.

    Owns the run's paths, shared by preparation, extraction and the pipeline.
    Construction does not run commands or create directories. run() checks
    the executable and working directory. Failure exceptions include only
    the last 100 lines per channel; the complete output remains in the log.
    """

    def __init__(self, working_dir, executable="gmx", logger=None):
        self.paths = RunPaths(working_dir)
        executable = str(executable)
        self.executable = str(Path(executable).expanduser().resolve()) if "/" in executable else executable
        self.logger = logger or logging.getLogger("g4_builder.gromacs")

    @property
    def working_dir(self):
        return self.paths.working_dir

    def check_available(self):
        if shutil.which(self.executable) is None:
            raise FileNotFoundError(f"GROMACS executable not found or not executable: {self.executable}")
        if not self.working_dir.is_dir():
            raise NotADirectoryError(self.working_dir)

    def run(self, arguments, input_text="", timeout=None):
        self.check_available()
        command = [self.executable, *map(str, arguments)]
        self.logger.debug("Running in %s: %s", self.working_dir, command)
        stdout, stderr = deque(maxlen=100), deque(maxlen=100)
        try:
            with subprocess.Popen(command, cwd=self.working_dir, text=True,
                                  encoding="utf-8", errors="replace", stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
                def consume(stream, name, tail):
                    for line in stream:
                        tail.append(line)
                        self.logger.debug("%s: %s", name, line.rstrip())
                threads = [threading.Thread(target=consume, args=(stream, name, tail), daemon=True)
                           for stream, name, tail in ((process.stdout, "stdout", stdout),
                                                       (process.stderr, "stderr", stderr))]
                for thread in threads:
                    thread.start()
                try:
                    try:
                        process.stdin.write(input_text)
                        process.stdin.close()
                    except BrokenPipeError:
                        pass
                    process.wait(timeout=timeout)
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                finally:
                    for thread in threads:
                        thread.join()
                if process.returncode:
                    raise subprocess.CalledProcessError(process.returncode, command,
                                                        "".join(stdout), "".join(stderr))
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
            self.logger.error("%s failed: %s", arguments[0], exc)
            for channel, tail in (("stdout", stdout), ("stderr", stderr)):
                if tail:
                    self.logger.error("%s:\n%s", channel, "".join(tail).rstrip())
            raise
        self.logger.info("%s completed successfully", arguments[0])

    def probe(self, arguments, timeout=2):
        """Capture a bounded interactive-menu probe; a timeout is expected."""
        self.check_available()
        command = [self.executable, *map(str, arguments)]
        stdbuf = shutil.which("stdbuf")
        if stdbuf:
            command = [stdbuf, "-oL", "-eL", *command]
        try:
            result = subprocess.run(command, cwd=self.working_dir, input="", text=True,
                                    capture_output=True, timeout=timeout, check=False)
            stdout, stderr = result.stdout, result.stderr
            code = result.returncode
        except subprocess.TimeoutExpired as exc:
            stdout, stderr = exc.stdout, exc.stderr
            code = 0
            self.logger.debug("Stopped menu probe after %s seconds", timeout)
        def decode(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else value or ""
        output = decode(stdout) + "\n" + decode(stderr)
        self.logger.debug("Force-field probe:\n%s", output)
        if code:
            self.logger.error("Force-field discovery failed:\n%s", output)
            raise subprocess.CalledProcessError(code, command, stdout, stderr)
        return output
