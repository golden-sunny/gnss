#!/usr/bin/env python3
"""Exercise 5.6 电离层全天绘图脚本。"""

from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots  # noqa: F401
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
SUMMARY_CSV = DATA_DIR / "exam-5.6-ionosphere_summary.csv"
DETAIL_CSV = DATA_DIR / "exam-5.6-ionosphere.csv"
OUTPUT_BASE = DATA_DIR / "exam-5.6-ionosphere_science"
G10_OUTPUT_BASE = DATA_DIR / "exam-5.6-ionosphere_g10_science"

warnings.filterwarnings(
    "ignore",
    message=r"Glyph .* missing from font\(s\) Times New Roman\.",
)


def get_english_fonts() -> list[str]:
    available = {font.name for font in font_manager.fontManager.ttflist}
    preferred = [
        "Times New Roman",
        "DejaVu Serif",
    ]
    return [font for font in preferred if font in available]


def get_chinese_font_name() -> str:
    available = {font.name for font in font_manager.fontManager.ttflist}
    preferred = [
        "SimSun",
        "STSong",
        "Noto Serif CJK SC",
        "Source Han Serif SC",
    ]
    selected = [font for font in preferred if font in available]
    if selected:
        return selected[0]
    # Fallback to a SciencePlots-compatible serif if no Chinese serif is available.
    return get_english_fonts()[0]


def configure_matplotlib() -> None:
    plt.style.use(["science", "ieee", "no-latex", "grid"])
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "font.family": "serif",
            "font.serif": get_english_fonts(),
            "axes.unicode_minus": False,
            "axes.linewidth": 1.0,
            "axes.labelsize": 9,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "grid.alpha": 0.45,
            "grid.linewidth": 0.6,
            "lines.linewidth": 1.4,
        }
    )


CN_FONT = FontProperties(family=get_chinese_font_name())
EN_FONT = FontProperties(family=get_english_fonts()[0])


def load_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["epoch"])
    numeric_cols = [
        "sat_count",
        "dual_count",
        "mean_elevation_deg",
        "mean_model_l1_iono_m",
        "mean_dual_l1_iono_m",
        "mean_diff_m",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values("epoch").reset_index(drop=True)


def load_detail(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["epoch"])
    numeric_cols = [
        "elevation_deg",
        "azimuth_deg",
        "leveled_dual_l1_iono_m",
        "klobuchar_l1_iono_m",
        "diff_m",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values(["sat", "epoch"]).reset_index(drop=True)


def insert_gap_rows(df: pd.DataFrame, time_col: str, value_cols: list[str]) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    df = df.sort_values(time_col).reset_index(drop=True)
    diffs = df[time_col].diff().dt.total_seconds()
    positive_diffs = diffs[diffs > 0]
    if positive_diffs.empty:
        return df

    nominal_step = positive_diffs.median()
    gap_threshold = max(nominal_step * 1.5, nominal_step + 1.0)

    rows: list[dict] = []
    for idx, row in df.iterrows():
        if idx > 0:
            dt = (row[time_col] - df.loc[idx - 1, time_col]).total_seconds()
            if dt > gap_threshold:
                gap_row = {col: pd.NA for col in df.columns}
                gap_row[time_col] = df.loc[idx - 1, time_col] + (row[time_col] - df.loc[idx - 1, time_col]) / 2
                for col in value_cols:
                    gap_row[col] = float("nan")
                rows.append(gap_row)
        rows.append(row.to_dict())

    return pd.DataFrame(rows, columns=df.columns)


def apply_time_axis(ax: plt.Axes, start, end) -> None:
    ax.set_xlim(start, end)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", linewidth=0.35, alpha=0.22)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(EN_FONT)


def save_figure(fig: plt.Figure, output_base: Path) -> tuple[Path, Path]:
    svg_path = output_base.parent / f"{output_base.name}.svg"
    png_path = output_base.parent / f"{output_base.name}.png"
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)
    return svg_path, png_path


def plot_g10(detail: pd.DataFrame) -> tuple[Path, Path]:
    g10 = detail.loc[detail["sat"] == "G10"].copy()
    if g10.empty:
        raise ValueError("未在 exam-5.6-ionosphere.csv 中找到 G10 卫星记录。")

    g10_plot = insert_gap_rows(
        g10,
        "epoch",
        ["klobuchar_l1_iono_m", "leveled_dual_l1_iono_m", "diff_m"],
    )

    start = g10["epoch"].min()
    end = g10["epoch"].max()
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(6.8, 4.8),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )

    ax_top.plot(
        g10_plot["epoch"],
        g10_plot["klobuchar_l1_iono_m"],
        color="black",
        linestyle="-",
        label="广播 Klobuchar 模型",
    )
    ax_top.plot(
        g10_plot["epoch"],
        g10_plot["leveled_dual_l1_iono_m"],
        color="#5a5a5a",
        linestyle="--",
        label="双频平滑反演值",
    )
    ax_top.set_ylabel("L1电离层延迟 / m", fontproperties=CN_FONT)
    ax_top.set_title("G10卫星电离层延迟", fontproperties=CN_FONT)
    ax_top.legend(loc="upper right", prop=CN_FONT)
    apply_time_axis(ax_top, start, end)

    ax_bottom.plot(
        g10_plot["epoch"],
        g10_plot["diff_m"],
        color="#303030",
        linestyle="-.",
        label="双频反演减模型值",
    )
    ax_bottom.axhline(0.0, color="#808080", linewidth=0.9)
    ax_bottom.set_ylabel("差值 / m", fontproperties=CN_FONT)
    ax_bottom.set_xlabel("GPST时间", fontproperties=CN_FONT)
    apply_time_axis(ax_bottom, start, end)

    fig.suptitle("练习5.6 G10卫星电离层延迟", x=0.08, y=1.02, ha="left", fontsize=10, fontweight="bold", fontproperties=CN_FONT)
    fig.text(
        0.08,
        0.985,
        (
            f"样本数 = {len(g10)}，高度角范围 = "
            f"{g10['elevation_deg'].min():.1f}°–{g10['elevation_deg'].max():.1f}°"
        ),
        ha="left",
        va="top",
        fontsize=8,
        fontproperties=CN_FONT,
    )
    return save_figure(fig, G10_OUTPUT_BASE)


def main() -> None:
    configure_matplotlib()
    summary = load_summary(SUMMARY_CSV)
    detail = load_detail(DETAIL_CSV)
    summary_plot = insert_gap_rows(
        summary,
        "epoch",
        ["mean_model_l1_iono_m", "mean_dual_l1_iono_m", "mean_diff_m"],
    )
    start = summary["epoch"].min()
    end = summary["epoch"].max()

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(6.8, 4.8),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )

    ax_top.plot(
        summary_plot["epoch"],
        summary_plot["mean_model_l1_iono_m"],
        color="black",
        linestyle="-",
        label="广播 Klobuchar 模型",
    )
    ax_top.plot(
        summary_plot["epoch"],
        summary_plot["mean_dual_l1_iono_m"],
        color="#5a5a5a",
        linestyle="--",
        label="双频平滑反演均值",
    )
    ax_top.set_ylabel("L1电离层延迟 / m", fontproperties=CN_FONT)
    ax_top.set_title("全天GPS L1电离层延迟", fontproperties=CN_FONT)
    ax_top.legend(loc="upper right", prop=CN_FONT)
    apply_time_axis(ax_top, start, end)

    ax_bottom.plot(
        summary_plot["epoch"],
        summary_plot["mean_diff_m"],
        color="#303030",
        linestyle="-.",
        label="双频反演减模型均值",
    )
    ax_bottom.axhline(0.0, color="#808080", linewidth=0.9)
    ax_bottom.set_ylabel("差值 / m", fontproperties=CN_FONT)
    ax_bottom.set_xlabel("GPST时间", fontproperties=CN_FONT)
    apply_time_axis(ax_bottom, start, end)

    fig.suptitle("练习5.6 电离层延迟全天变化", x=0.08, y=1.02, ha="left", fontsize=10, fontweight="bold", fontproperties=CN_FONT)
    fig.text(
        0.08,
        0.985,
        (
            f"历元数 = {len(summary)}，平均可见GPS卫星数 = {summary['sat_count'].mean():.1f}，"
            f"平均双频有效卫星数 = {summary['dual_count'].mean():.1f}"
        ),
        ha="left",
        va="top",
        fontsize=8,
        fontproperties=CN_FONT,
    )

    svg_path, png_path = save_figure(fig, OUTPUT_BASE)
    plt.close(fig)
    g10_svg_path, g10_png_path = plot_g10(detail)
    plt.close("all")

    print(f"SciencePlots SVG written to: {svg_path}")
    print(f"SciencePlots PNG written to: {png_path}")
    print(f"G10 SciencePlots SVG written to: {g10_svg_path}")
    print(f"G10 SciencePlots PNG written to: {g10_png_path}")


if __name__ == "__main__":
    main()
