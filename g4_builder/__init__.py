"""DNA folding preparation library. Importing never starts GROMACS."""

import logging as _logging

_logging.getLogger(__name__).addHandler(_logging.NullHandler())

from .config import SimulationConfig, ScheduleConfig
from .paths import RunPaths
from .pipeline import FoldingPipeline, RunResult
from .strand import DNAStrandBuilder
from .references import build_references, Reference, ReferenceSet
from .template_builder import build_reference_template
from .protocol import Protocol, build_protocol
from .runner import GromacsRunner
from .gromacs import GromacsPreparation, SimulationOutputs
from .extraction import TrajectoryFrameExtractor
from .logging import configure_logging
from .pdb_io import read_pdb, save_pdb

__all__ = [
    "SimulationConfig",
    "ScheduleConfig",
    "RunPaths",
    "FoldingPipeline",
    "RunResult",
    "DNAStrandBuilder",
    "build_references",
    "build_reference_template",
    "Reference",
    "ReferenceSet",
    "Protocol",
    "build_protocol",
    "GromacsRunner",
    "GromacsPreparation",
    "SimulationOutputs",
    "TrajectoryFrameExtractor",
    "configure_logging",
    "read_pdb",
    "save_pdb",
]
