"""BioPilot - A lightweight, reproducible bioinformatics assistant."""

__version__ = "0.1.0"
__author__ = "BioPilot"

from biopilot.src.fetcher import DatasetFetcher
from biopilot.src.annotation import AnnotationDB
from biopilot.src.pipeline import PipelineManager
from biopilot.src.analyzer import DataAnalyzer
from biopilot.src.reproducibility import ReproducibilityLogger
