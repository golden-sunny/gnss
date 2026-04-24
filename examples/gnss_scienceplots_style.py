#!/usr/bin/env python3
"""Shared plotting helpers for GNSS paper-style SciencePlots figures."""

from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots  # noqa: F401
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties


warnings.filterwarnings(
    "ignore",
    message=r"Glyph .* missing from font\(s\) Times New Roman\.",
)


def _available_font_names() -> set[str]:
    return {font.name for font in font_manager.fontManager.ttflist}


def get_english_font_name() -> str:
    available = _available_font_names()
    for name in ["Times New Roman", "DejaVu Serif"]:
        if name in available:
            return name
    return "DejaVu Serif"


def get_chinese_font_name() -> str:
    available = _available_font_names()
    for name in ["SimSun", "STSong", "Noto Serif CJK SC", "Source Han Serif SC"]:
        if name in available:
            return name
    return get_english_font_name()


CN_FONT = FontProperties(family=get_chinese_font_name())
EN_FONT = FontProperties(family=get_english_font_name())


def configure_scienceplots(compact: bool = True) -> None:
    styles = ["science", "ieee", "no-latex"] if compact else ["science", "no-latex", "grid"]
    plt.style.use(styles)
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "font.family": "serif",
            "font.serif": [get_english_font_name()],
            "axes.unicode_minus": False,
            "axes.linewidth": 1.0,
            "axes.facecolor": "white",
            "axes.edgecolor": "#444444",
            "axes.labelsize": 9 if compact else 12,
            "axes.titlesize": 9.5 if compact else 14,
            "axes.titleweight": "bold",
            "xtick.labelsize": 8 if compact else 10.5,
            "ytick.labelsize": 8 if compact else 10.5,
            "legend.fontsize": 8 if compact else 10.5,
            "grid.color": "#d7d7d7",
            "grid.linewidth": 0.6 if compact else 0.7,
            "grid.alpha": 0.6 if compact else 0.8,
            "lines.linewidth": 1.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "legend.frameon": True,
            "legend.facecolor": "white",
            "legend.edgecolor": "#c8c8c8",
            "legend.framealpha": 1.0,
        }
    )


def apply_common_axis_style(ax: plt.Axes) -> None:
    ax.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.7)
    ax.grid(True, which="minor", linestyle="--", linewidth=0.35, alpha=0.28)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_color("#444444")
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(EN_FONT)


def apply_time_axis(
    ax: plt.Axes,
    start,
    end,
    hour_interval: int = 3,
    hour_only: bool = False,
) -> None:
    ax.set_xlim(start, end)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
    if hour_only:
        ax.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda value, _: str(mdates.num2date(value).hour))
        )
    else:
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    apply_common_axis_style(ax)


def save_figure(fig: plt.Figure, output_base: Path) -> tuple[Path, Path]:
    svg_path = output_base.parent / f"{output_base.name}.svg"
    png_path = output_base.parent / f"{output_base.name}.png"
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)
    return svg_path, png_path


def insert_gap_rows(df: pd.DataFrame, time_col: str, value_cols: list[str]) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    ordered = df.sort_values(time_col).reset_index(drop=True)
    diffs = ordered[time_col].diff().dt.total_seconds()
    positive_diffs = diffs[diffs > 0]
    if positive_diffs.empty:
        return ordered

    nominal_step = positive_diffs.median()
    gap_threshold = max(nominal_step * 1.5, nominal_step + 1.0)

    rows: list[dict] = []
    for idx, row in ordered.iterrows():
        if idx > 0:
            dt = (row[time_col] - ordered.loc[idx - 1, time_col]).total_seconds()
            if dt > gap_threshold:
                gap_row = {col: pd.NA for col in ordered.columns}
                gap_row[time_col] = ordered.loc[idx - 1, time_col] + (
                    row[time_col] - ordered.loc[idx - 1, time_col]
                ) / 2
                for col in value_cols:
                    gap_row[col] = float("nan")
                rows.append(gap_row)
        rows.append(row.to_dict())

    return pd.DataFrame(rows, columns=ordered.columns)


def muted_palette() -> list[str]:
    return [
        "#0C5DA5",
        "#00B945",
        "#FF9500",
        "#845B97",
        "#474747",
        "#9E9E9E",
        "#5C7CFA",
        "#2B8A3E",
        "#C77D00",
        "#7A5C99",
    ]
