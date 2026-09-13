"""Shared manuscript figure styling for the manuscript package."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl


# A restrained editorial-print palette. Colour separates biological categories,
# while neutral paper-like grays, square annotations, and fine rules avoid a
# dashboard or presentation-template appearance.
INK = "#1F2020"
MUTED = "#626461"
GRID = "#E8E6E1"
BLUE = "#506573"
SKY = "#829EA3"
GREEN = "#526E60"
ORANGE = "#98664F"
RED = "#8D4D52"
PURPLE = "#6C5D70"
GOLD = "#99834D"
GRAY = "#858783"
LIGHT = "#EFEDE8"


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.4,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.0,
            "axes.titleweight": "normal",
            "axes.edgecolor": MUTED,
            "axes.linewidth": 0.7,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.labelsize": 7.7,
            "ytick.labelsize": 7.7,
            "legend.fontsize": 7.4,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_title(ax, letter: str, title: str) -> None:
    ax.set_title(title, loc="left", pad=6)
    ax.text(
        -0.12,
        1.035,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        color=INK,
        va="bottom",
    )


def finish_axis(ax, *, grid: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.55, alpha=0.75)
        ax.set_axisbelow(True)


def save_figure(fig, png: Path, pdf: Path | None = None) -> None:
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=450)
    fig.savefig(pdf or png.with_suffix(".pdf"))
