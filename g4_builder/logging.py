"""Opt-in logging configuration for the library's calling application."""

import logging
from pathlib import Path

LOGGER_NAME = "g4_builder"


class _OutputFormatter(logging.Formatter):
    """Keep subprocess text intact and format only application messages."""

    def __init__(self, fmt):
        super().__init__(fmt)
        self._at_line_start = True

    def format(self, record):
        if getattr(record, "gromacs_output", False):
            text = record.getMessage()
        else:
            text = super().format(record) + "\n"
            # Separate status messages from output without a final newline.
            if not self._at_line_start:
                text = "\n" + text
        if text:
            self._at_line_start = text.endswith(("\n", "\r"))
        return text


def configure_logging(log_file=None, console=True, console_level=logging.INFO):
    """Configure package-owned handlers; repeated calls replace them safely.

    Command output is unprefixed DEBUG text in the file, under a block header.
    The console displays brief INFO status and errors. Without this function,
    normal Python logging
    configuration controls the library; importing it creates no files.
    """
    logger = logging.getLogger(LOGGER_NAME)
    # Leave application-owned handlers alone when replacing our configuration.
    for handler in list(logger.handlers):
        if getattr(handler, "_g4_owned", False):
            logger.removeHandler(handler)
            handler.close()
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    handlers = []
    if console:
        handler = logging.StreamHandler()
        handler.setLevel(console_level)
        handler.setFormatter(_OutputFormatter("%(levelname)s: %(message)s"))
        handlers.append(handler)
    if log_file is not None:
        path = Path(log_file).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(path, encoding="utf-8")
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(
            _OutputFormatter("%(asctime)s %(levelname)s %(message)s")
        )
        handlers.append(handler)
    for handler in handlers:
        # The formatter supplies newlines for messages; raw output has its own.
        handler.terminator = ""
        handler._g4_owned = True
        logger.addHandler(handler)
    return logger
