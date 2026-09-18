from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import pandas as pd

from analytics.tracking import (
    ANALYZE_DATASET,
    SUMMARY_DATASET,
    FIGURES_DIR,
    rebuild_category_dataset,
    rebuild_summary_dataset,
)

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
LOG = logging.getLogger(__name__)


def _load_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        LOG.warning("Dataset not found: %s", path)
        return None
    return pd.read_csv(path)


def _plot_bar(series, title: str, xlabel: str, destination: Path) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ax = series.plot(kind="bar", figsize=(8, 5))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(destination)
    plt.close()
    LOG.info("Wrote %s", destination)


def _plot_horizontal(series, title: str, xlabel: str, destination: Path) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ax = series.sort_values().plot(kind="barh", figsize=(8, 6))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    plt.tight_layout()
    plt.savefig(destination)
    plt.close()
    LOG.info("Wrote %s", destination)


def build_summary_charts(summary_df: pd.DataFrame) -> None:
    if summary_df.empty:
        LOG.warning("Summary dataset is empty; skipping charts.")
        return
    _plot_bar(
        summary_df["dominant_pathogen_type"].value_counts(dropna=False),
        title="Dominant pathogen type counts",
        xlabel="Pathogen type",
        destination=FIGURES_DIR / "dominant_pathogen_type_counts.png",
    )
    _plot_horizontal(
        summary_df.set_index("patient_id")["overall_confidence"],
        title="Overall confidence by patient",
        xlabel="Confidence",
        destination=FIGURES_DIR / "overall_confidence_by_patient.png",
    )


def build_category_charts(category_df: pd.DataFrame) -> None:
    if category_df.empty:
        LOG.warning("Category dataset is empty; skipping charts.")
        return
    _plot_bar(
        category_df["category"].value_counts(dropna=False),
        title="Agent category counts",
        xlabel="Category",
        destination=FIGURES_DIR / "agent_category_counts.png",
    )


def rebuild_datasets(inputs: Sequence[Path]) -> None:
    rebuild_category_dataset(inputs)
    rebuild_summary_dataset(inputs)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate analytics charts from CSV datasets.")
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Rebuild analytics CSV datasets by scanning patient folders.",
    )
    parser.add_argument(
        "--inputs",
        nargs="*",
        type=Path,
        help="Patient folders or parent directories (required when using --rebuild).",
    )
    parser.add_argument(
        "--skip-summary",
        action="store_true",
        help="Skip summary charts.",
    )
    parser.add_argument(
        "--skip-category",
        action="store_true",
        help="Skip category charts.",
    )
    args = parser.parse_args()

    if args.rebuild:
        if not args.inputs:
            parser.error("--rebuild requires at least one path via --inputs")
        rebuild_datasets(args.inputs)

    summary_df = None if args.skip_summary else _load_csv(SUMMARY_DATASET)
    category_df = None if args.skip_category else _load_csv(ANALYZE_DATASET)

    if summary_df is not None and not args.skip_summary:
        build_summary_charts(summary_df)
    if category_df is not None and not args.skip_category:
        build_category_charts(category_df)


if __name__ == "__main__":
    main()

