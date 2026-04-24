#!/usr/bin/env python3
"""Unified GNSS SciencePlots figures for Exercise 5.8."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from gnss_scienceplots_style import (
    CN_FONT,
    apply_time_axis,
    configure_scienceplots,
    insert_gap_rows,
    save_figure,
)


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
SUMMARY_CSV = DATA_DIR / "exam-5.8-tropo_iono_satpos_summary.csv"
DETAIL_CSV = DATA_DIR / "exam-5.8-tropo_iono_satpos.csv"

OVERVIEW_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos"
G01_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos_g01"
OVERVIEW_IONO_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos_iono_science_ieee"
OVERVIEW_TROPO_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos_tropo_science_ieee"
G01_IONO_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos_g01_iono_science_ieee"
G01_TROPO_BASE = DATA_DIR / "exam-5.8-tropo_iono_satpos_g01_tropo_science_ieee"


def load_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["epoch"])
    for col in [
        "mean_iono_model_m",
        "mean_dual_iono_m",
        "mean_trop_m",
        "mean_elevation_deg",
        "sat_count",
        "dual_count",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values("epoch").reset_index(drop=True)


def load_detail(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["epoch"])
    for col in [
        "model_iono_m",
        "model_trop_m",
        "dual_iono_m",
        "dual_minus_model_m",
        "elevation_deg",
        "azimuth_deg",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values(["sat", "epoch"]).reset_index(drop=True)


def finalize_overview(fig: plt.Figure, title: str, subtitle: str) -> None:
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=16, fontweight="bold", fontproperties=CN_FONT)
    fig.text(0.08, 0.95, subtitle, ha="left", va="top", fontsize=10.5, color="#4f4f4f", fontproperties=CN_FONT)


def plot_overview(summary: pd.DataFrame, output_base: Path) -> tuple[Path, Path]:
    summary_plot = insert_gap_rows(
        summary,
        "epoch",
        ["mean_iono_model_m", "mean_dual_iono_m", "mean_trop_m"],
    )
    start = summary["epoch"].min()
    end = summary["epoch"].max()

    configure_scienceplots(compact=False)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.8, 8.2), sharex=True, constrained_layout=False)

    ax1.plot(summary_plot["epoch"], summary_plot["mean_iono_model_m"], color="#0C5DA5", label="广播模型电离层延迟")
    ax1.plot(summary_plot["epoch"], summary_plot["mean_dual_iono_m"], color="#00B945", linestyle="--", label="双频平滑电离层延迟")
    ax1.set_ylabel("延迟 / m", fontproperties=CN_FONT)
    ax1.set_title("历元平均电离层延迟", fontproperties=CN_FONT)
    ax1.legend(loc="upper right", prop=CN_FONT)
    ax1.text(0.0, 1.03, "(a)", transform=ax1.transAxes, fontsize=12.5, fontweight="bold")
    apply_time_axis(ax1, start, end)
    ax1.tick_params(axis="x", labelbottom=True)

    ax2.plot(summary_plot["epoch"], summary_plot["mean_trop_m"], color="#222222", linestyle="-.", label="Saastamoinen 对流层延迟")
    ax2.set_ylabel("延迟 / m", fontproperties=CN_FONT)
    ax2.set_xlabel("时间 / h", fontproperties=CN_FONT)
    ax2.set_title("历元平均对流层延迟", fontproperties=CN_FONT)
    ax2.legend(loc="upper right", prop=CN_FONT)
    ax2.text(0.0, 1.03, "(b)", transform=ax2.transAxes, fontsize=12.5, fontweight="bold")
    apply_time_axis(ax2, start, end)

    finalize_overview(
        fig,
        "练习5.8 电离层与对流层延迟总体变化",
        f"历元数 = {len(summary)}，平均可见GPS卫星数 = {summary['sat_count'].mean():.1f}，平均双频有效卫星数 = {summary['dual_count'].mean():.1f}",
    )
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.10, hspace=0.30)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def plot_g01(detail: pd.DataFrame, output_base: Path) -> tuple[Path, Path]:
    g01 = detail.loc[detail["sat"] == "G01"].copy()
    if g01.empty:
        raise ValueError("未在明细 CSV 中找到 G01 卫星记录。")

    g01_plot = insert_gap_rows(g01, "epoch", ["model_iono_m", "model_trop_m", "dual_iono_m"])
    start = g01["epoch"].min()
    end = g01["epoch"].max()

    configure_scienceplots(compact=False)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.8, 8.2), sharex=True, constrained_layout=False)

    ax1.plot(g01_plot["epoch"], g01_plot["model_iono_m"], color="#0C5DA5", label="广播模型电离层延迟")
    ax1.plot(g01_plot["epoch"], g01_plot["dual_iono_m"], color="#00B945", linestyle="--", label="双频平滑电离层延迟")
    ax1.set_ylabel("延迟 / m", fontproperties=CN_FONT)
    ax1.set_title("G01卫星电离层延迟", fontproperties=CN_FONT)
    ax1.legend(loc="upper right", prop=CN_FONT)
    ax1.text(0.0, 1.03, "(a)", transform=ax1.transAxes, fontsize=12.5, fontweight="bold")
    apply_time_axis(ax1, start, end)
    ax1.tick_params(axis="x", labelbottom=True)

    ax2.plot(g01_plot["epoch"], g01_plot["model_trop_m"], color="#222222", linestyle="-.", label="Saastamoinen 对流层延迟")
    ax2.set_ylabel("延迟 / m", fontproperties=CN_FONT)
    ax2.set_xlabel("时间 / h", fontproperties=CN_FONT)
    ax2.set_title("G01卫星对流层延迟", fontproperties=CN_FONT)
    ax2.legend(loc="upper right", prop=CN_FONT)
    ax2.text(0.0, 1.03, "(b)", transform=ax2.transAxes, fontsize=12.5, fontweight="bold")
    apply_time_axis(ax2, start, end)

    valid_dual = g01["dual_iono_m"].notna().sum()
    subtitle = (
        f"样本数 = {len(g01)}，双频有效样本数 = {valid_dual}，"
        f"高度角范围 = {g01['elevation_deg'].min():.1f}°–{g01['elevation_deg'].max():.1f}°"
    )
    finalize_overview(fig, "练习5.8 G01卫星延迟变化", subtitle)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.10, hspace=0.30)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def plot_single_panel(
    x_data: pd.Series,
    y_series: list[tuple[pd.Series, str, str, str]],
    output_base: Path,
    y_label: str,
    x_label: str,
    start,
    end,
) -> tuple[Path, Path]:
    configure_scienceplots(compact=True)
    fig, ax = plt.subplots(figsize=(3.3, 2.5), constrained_layout=True)
    for y_data, label, color, linestyle in y_series:
        ax.plot(x_data, y_data, color=color, linestyle=linestyle, label=label)
    ax.set_ylabel(y_label, fontproperties=CN_FONT)
    ax.set_xlabel(x_label, fontproperties=CN_FONT)
    ax.legend(loc="best", prop=CN_FONT)
    apply_time_axis(ax, start, end, hour_interval=3, hour_only=True)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def main() -> None:
    summary = load_summary(SUMMARY_CSV)
    detail = load_detail(DETAIL_CSV)

    overview_svg, overview_png = plot_overview(summary, OVERVIEW_BASE)
    g01_svg, g01_png = plot_g01(detail, G01_BASE)

    summary_plot = insert_gap_rows(summary, "epoch", ["mean_iono_model_m", "mean_dual_iono_m", "mean_trop_m"])
    g01 = detail.loc[detail["sat"] == "G01"].copy()
    g01_plot = insert_gap_rows(g01, "epoch", ["model_iono_m", "dual_iono_m", "model_trop_m"])
    start = summary["epoch"].min()
    end = summary["epoch"].max()

    overview_iono = plot_single_panel(
        summary_plot["epoch"],
        [
            (summary_plot["mean_iono_model_m"], "广播模型电离层延迟", "#0C5DA5", "-"),
            (summary_plot["mean_dual_iono_m"], "双频平滑电离层延迟", "#00B945", "--"),
        ],
        OVERVIEW_IONO_BASE,
        "延迟 / m",
        "时间 / h",
        start,
        end,
    )
    overview_tropo = plot_single_panel(
        summary_plot["epoch"],
        [(summary_plot["mean_trop_m"], "Saastamoinen 对流层延迟", "#222222", "-.")],
        OVERVIEW_TROPO_BASE,
        "延迟 / m",
        "时间 / h",
        start,
        end,
    )
    g01_iono = plot_single_panel(
        g01_plot["epoch"],
        [
            (g01_plot["model_iono_m"], "广播模型电离层延迟", "#0C5DA5", "-"),
            (g01_plot["dual_iono_m"], "双频平滑电离层延迟", "#00B945", "--"),
        ],
        G01_IONO_BASE,
        "延迟 / m",
        "时间 / h",
        g01["epoch"].min(),
        g01["epoch"].max(),
    )
    g01_tropo = plot_single_panel(
        g01_plot["epoch"],
        [(g01_plot["model_trop_m"], "Saastamoinen 对流层延迟", "#222222", "-.")],
        G01_TROPO_BASE,
        "延迟 / m",
        "时间 / h",
        g01["epoch"].min(),
        g01["epoch"].max(),
    )

    print(f"总体图 SVG 已输出: {overview_svg}")
    print(f"总体图 PNG 已输出: {overview_png}")
    print(f"G01 图 SVG 已输出: {g01_svg}")
    print(f"G01 图 PNG 已输出: {g01_png}")
    for label, outputs in [
        ("SciencePlots 总体电离层图", overview_iono),
        ("SciencePlots 总体对流层图", overview_tropo),
        ("SciencePlots G01 电离层图", g01_iono),
        ("SciencePlots G01 对流层图", g01_tropo),
    ]:
        print(f"{label} SVG 已输出: {outputs[0]}")
        print(f"{label} PNG 已输出: {outputs[1]}")


if __name__ == "__main__":
    main()
