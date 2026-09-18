"""RAG_re: a small, replayable PubMed ablation pipeline."""

from .engine import RagReEngine
from .versioning import PIPELINE_VERSION

__all__ = ["RagReEngine", "PIPELINE_VERSION"]
__version__ = PIPELINE_VERSION
