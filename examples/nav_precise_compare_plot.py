#!/usr/bin/env python3
"""Unified GNSS SciencePlots figures for nav_precise_compare outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from gnss_scienceplots_style import (
    CN_FONT,
    EN_FONT,
    apply_time_axis,
    configure_scienceplots,
    muted_palette,
    save_figure,
)


ROOT = Path(__file__).resolve().parents[1]
DAILY_CSV = ROOT / "nav_precise_compare_20250101_30s_daily.csv"
EPOCH_CSV = ROOT / "nav_precise_compare_20250101_000500.csv"


def load_daily(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (TypeError, ValueError):
            pass
    df["epoch"] = pd.Timestamp("2025-01-01") + pd.to_timedelta(df["sec_of_day"], unit="s")
    return df.sort_values("epoch").reset_index(drop=True)


def load_epoch(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ["pos_dx", "pos_dy", "pos_dz", "vel_dx", "vel_dy", "vel_dz", "pos_norm", "vel_norm"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values(["system", "sat"]).reset_index(drop=True)


def plot_daily_components(
    daily: pd.DataFrame,
    system: str,
    prefix: str,
    title: str,
    y_label: str,
    output_base: Path,
) -> tuple[Path, Path]:
    configure_scienceplots(compact=False)
    start = daily["epoch"].min()
    end = daily["epoch"].max()
    fig, axes = plt.subplots(3, 1, figsize=(10.6, 7.8), sharex=True, constrained_layout=False)
    colors = ["#0C5DA5", "#00B945", "#222222"]
    component_titles = ["X分量", "Y分量", "Z分量"]

    for ax, component, color, component_title in zip(axes, ["x", "y", "z"], colors, component_titles):
        ax.plot(daily["epoch"], daily[f"{prefix}_{component}"], color=color, linewidth=1.8)
        ax.set_ylabel(y_label, fontproperties=CN_FONT)
        ax.set_title(component_title, loc="left", fontproperties=CN_FONT)
        apply_time_axis(ax, start, end)

    axes[-1].set_xlabel("时间 / h", fontproperties=CN_FONT)
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=15, fontweight="bold", fontproperties=CN_FONT)
    fig.text(
        0.08,
        0.95,
        f"{system} 系统；基于 2025-01-01 全天 30 s 统计结果。",
        ha="left",
        va="top",
        fontsize=10,
        color="#4f4f4f",
        fontproperties=CN_FONT,
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.88, bottom=0.10, hspace=0.28)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def plot_epoch_norms(epoch_df: pd.DataFrame, system: str, output_base: Path) -> tuple[Path, Path]:
    configure_scienceplots(compact=False)
    subset = epoch_df.loc[epoch_df["system"] == system].copy()
    if subset.empty:
        raise ValueError(f"未找到 {system} 系统的单历元数据。")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5, 6.8), sharex=True, constrained_layout=False)
    x = range(len(subset))
    colors = muted_palette()
    bar_colors = [colors[idx % len(colors)] for idx in x]

    ax1.bar(x, subset["pos_norm"], color=bar_colors, edgecolor="#444444", linewidth=0.4)
    ax1.set_ylabel("位置范数 / m", fontproperties=CN_FONT)
    ax1.set_title("单历元位置差范数", loc="left", fontproperties=CN_FONT)
    ax1.grid(True, axis="y", linestyle="--", linewidth=0.7, alpha=0.7)

    ax2.bar(x, subset["vel_norm"], color=bar_colors, edgecolor="#444444", linewidth=0.4)
    ax2.set_ylabel("速度范数 / m/s", fontproperties=CN_FONT)
    ax2.set_title("单历元速度差范数", loc="left", fontproperties=CN_FONT)
    ax2.set_xlabel("卫星 PRN", fontproperties=CN_FONT)
    ax2.grid(True, axis="y", linestyle="--", linewidth=0.7, alpha=0.7)

    ax2.set_xticks(list(x))
    ax2.set_xticklabels(subset["sat"], rotation=90)
    for ax in (ax1, ax2):
        for label in ax.get_xticklabels():
            label.set_fontproperties(EN_FONT)
        for label in ax.get_yticklabels():
            label.set_fontproperties(EN_FONT)
        ax.minorticks_on()

    fig.suptitle(
        f"{system} 广播星历与精密星历单历元差值范数",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=15,
        fontweight="bold",
        fontproperties=CN_FONT,
    )
    fig.text(
        0.08,
        0.95,
        f"历元 = 2025-01-01 00:05:00；卫星数量 = {len(subset)}。",
        ha="left",
        va="top",
        fontsize=10,
        color="#4f4f4f",
        fontproperties=CN_FONT,
    )
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.18, hspace=0.30)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def main() -> None:
    daily = load_daily(DAILY_CSV)
    epoch_df = load_epoch(EPOCH_CSV)

    outputs = [
        ("GPS日变化位置差图", plot_daily_components(
            daily, "GPS", "gps_pos", "GPS广播星历与精密星历位置差日变化", "位置差 / m", ROOT / "nav_precise_compare_20250101_30s_gps_pos_diff_science"
        )),
        ("GPS日变化速度差图", plot_daily_components(
            daily, "GPS", "gps_vel", "GPS广播星历与精密星历速度差日变化", "速度差 / m/s", ROOT / "nav_precise_compare_20250101_30s_gps_vel_diff_science"
        )),
        ("北斗日变化位置差图", plot_daily_components(
            daily, "BDS", "bds_pos", "北斗广播星历与精密星历位置差日变化", "位置差 / m", ROOT / "nav_precise_compare_20250101_30s_bds_pos_diff_science"
        )),
        ("北斗日变化速度差图", plot_daily_components(
            daily, "BDS", "bds_vel", "北斗广播星历与精密星历速度差日变化", "速度差 / m/s", ROOT / "nav_precise_compare_20250101_30s_bds_vel_diff_science"
        )),
        ("GPS单历元范数图", plot_epoch_norms(epoch_df, "GPS", ROOT / "nav_precise_compare_20250101_000500_gps_norm_science")),
        ("北斗单历元范数图", plot_epoch_norms(epoch_df, "BDS", ROOT / "nav_precise_compare_20250101_000500_bds_norm_science")),
    ]

    for label, (svg_path, png_path) in outputs:
        print(f"{label} SVG 已输出: {svg_path}")
        print(f"{label} PNG 已输出: {png_path}")


if __name__ == "__main__":
    main()
