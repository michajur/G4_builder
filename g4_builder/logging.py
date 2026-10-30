"""Opt-in logging configuration for the library's calling application."""
import logging
from pathlib import Path

LOGGER_NAME = "g4_builder"


def configure_logging(log_file=None, console=True, console_level=logging.INFO):
    """Configure package-owned handlers; repeated calls replace them safely.

    All command output is DEBUG in the file. The console displays brief INFO
    status and errors. With no call to this function, normal Python logging
    configuration controls the library; importing it creates no files.
    """
    logger = logging.getLogger(LOGGER_NAME)
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
        handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        handlers.append(handler)
    if log_file is not None:
        path = Path(log_file).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(path, encoding="utf-8")
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        handlers.append(handler)
    for handler in handlers:
        handler._g4_owned = True
        logger.addHandler(handler)
    return logger
