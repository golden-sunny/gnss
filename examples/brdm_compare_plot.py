#!/usr/bin/env python3
"""Unified GNSS SciencePlots figures for brdm_compare outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from gnss_scienceplots_style import (
    CN_FONT,
    apply_time_axis,
    configure_scienceplots,
    insert_gap_rows,
    muted_palette,
    save_figure,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "result" / "nav_sp3_compare"
CSV_FILE = RESULT_DIR / "brdm_compare_20250101_30s_satellite_diff.csv"


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ["sec_of_day", "pos_dx", "pos_dy", "pos_dz", "vel_dx", "vel_dy", "vel_dz"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["epoch"] = pd.Timestamp("2025-01-01") + pd.to_timedelta(df["sec_of_day"], unit="s")
    return df.sort_values(["system", "satellite", "epoch"]).reset_index(drop=True)


def plot_component_figure(
    data: pd.DataFrame,
    system: str,
    value_cols: list[str],
    title: str,
    y_label: str,
    output_base: Path,
) -> tuple[Path, Path]:
    configure_scienceplots(compact=False)
    subset = data.loc[data["system"] == system].copy()
    if subset.empty:
        raise ValueError(f"未找到 {system} 系统的数据。")

    start = subset["epoch"].min()
    end = subset["epoch"].max()
    satellites = sorted(subset["satellite"].unique())
    palette = muted_palette()

    fig, axes = plt.subplots(3, 1, figsize=(11.5, 8.0), sharex=True, constrained_layout=False)
    component_titles = ["X分量", "Y分量", "Z分量"]

    for sat_idx, satellite in enumerate(satellites):
        sat_df = subset.loc[subset["satellite"] == satellite].copy()
        sat_plot = insert_gap_rows(sat_df, "epoch", value_cols)
        linestyle = "--" if (system == "BDS" and sat_df["class"].iloc[0] == "GEO") else "-"
        color = palette[sat_idx % len(palette)]
        for ax, component_title, col in zip(axes, component_titles, value_cols):
            ax.plot(
                sat_plot["epoch"],
                sat_plot[col],
                color=color,
                linestyle=linestyle,
                linewidth=1.0,
                label=satellite,
            )
            ax.set_ylabel(y_label, fontproperties=CN_FONT)
            ax.set_title(component_title, loc="left", fontproperties=CN_FONT)
            apply_time_axis(ax, start, end)

    axes[-1].set_xlabel("时间 / h", fontproperties=CN_FONT)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.90, 0.5),
        ncol=2,
        prop=CN_FONT,
        frameon=True,
    )
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=15, fontweight="bold", fontproperties=CN_FONT)
    subtitle = f"{system} 卫星数量 = {len(satellites)}。"
    if system == "BDS":
        subtitle += " GEO 轨道使用虚线表示。"
    fig.text(
        0.08,
        0.95,
        subtitle,
        ha="left",
        va="top",
        fontsize=10,
        color="#4f4f4f",
        fontproperties=CN_FONT,
    )
    fig.subplots_adjust(left=0.10, right=0.86, top=0.88, bottom=0.10, hspace=0.28)
    outputs = save_figure(fig, output_base)
    plt.close(fig)
    return outputs


def main() -> None:
    data = load_data(CSV_FILE)

    outputs = [
        ("GPS位置差图", plot_component_figure(
            data,
            "GPS",
            ["pos_dx", "pos_dy", "pos_dz"],
            "GPS广播星历与精密星历位置差",
            "位置差 / m",
            RESULT_DIR / "brdm_compare_20250101_30s_gps_pos_diff_science",
        )),
        ("GPS速度差图", plot_component_figure(
            data,
            "GPS",
            ["vel_dx", "vel_dy", "vel_dz"],
            "GPS广播星历与精密星历速度差",
            "速度差 / m/s",
            RESULT_DIR / "brdm_compare_20250101_30s_gps_vel_diff_science",
        )),
        ("北斗位置差图", plot_component_figure(
            data,
            "BDS",
            ["pos_dx", "pos_dy", "pos_dz"],
            "北斗广播星历与精密星历位置差",
            "位置差 / m",
            RESULT_DIR / "brdm_compare_20250101_30s_bds_pos_diff_science",
        )),
        ("北斗速度差图", plot_component_figure(
            data,
            "BDS",
            ["vel_dx", "vel_dy", "vel_dz"],
            "北斗广播星历与精密星历速度差",
            "速度差 / m/s",
            RESULT_DIR / "brdm_compare_20250101_30s_bds_vel_diff_science",
        )),
    ]

    for label, (svg_path, png_path) in outputs:
        print(f"{label} SVG 已输出: {svg_path}")
        print(f"{label} PNG 已输出: {png_path}")


if __name__ == "__main__":
    main()
