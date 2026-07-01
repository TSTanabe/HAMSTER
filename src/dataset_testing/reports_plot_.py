#!/usr/bin/env python3

import argparse
import ast
import os
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


LABELS = ["TP", "FP", "FN", "TN"]
COLORS = {
    "TP": "#4daf4a",
    "FP": "#e41a1c",
    "FN": "#ff7f00",
    "TN": "#999999",
}


def parse_neighborhood(value):
    if pd.isna(value):
        return "singleton"

    if isinstance(value, str):
        try:
            value = ast.literal_eval(value)
        except Exception:
            return value

    if not value:
        return "singleton"

    genes = []
    for domains in value:
        if not domains:
            genes.append("no_neighbours")
        elif isinstance(domains, (list, tuple)):
            genes.append("-".join(str(x) for x in domains))
        else:
            genes.append(str(domains))

    return "/".join(genes)


def normalize_above_noise(value):
    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return True

    value = str(value).strip().lower()

    if value in {"true", "1", "yes", "above_noise", "above"}:
        return True

    if value in {"false", "0", "no", "below_noise", "below"}:
        return False

    return True


def classify_row(row, cutoff):
    # Below-noise extension rows are always counted as TN.
    if "above_noise" in row and not normalize_above_noise(row["above_noise"]):
        return "TN"

    if "below_noise_extension" in row:
        val = row["below_noise_extension"]
        if isinstance(val, bool) and val:
            return "TN"
        if str(val).strip().lower() in {"true", "1", "yes"}:
            return "TN"

    true_value = str(row["training_assignment"]).strip()
    predicted_positive = float(row["bitscore"]) >= float(cutoff)

    if predicted_positive and true_value == "TP":
        return "TP"
    if predicted_positive and true_value == "TN":
        return "FP"
    if (not predicted_positive) and true_value == "TP":
        return "FN"
    return "TN"


def load_cutoffs(cutoff_file, protein_name):
    df = pd.read_csv(cutoff_file, sep="\t")

    name_col = "name"
    if name_col not in df.columns:
        raise ValueError(f"Cutoff file needs a '{name_col}' column")

    row = df[df[name_col].astype(str) == protein_name]

    if row.empty:
        # fallback: detailed file may be grp2_DsrL_detailed.txt, cutoff name may be grp2_DsrL
        short_name = protein_name.replace("_detailed", "")
        row = df[df[name_col].astype(str) == short_name]

    if row.empty:
        raise ValueError(f"No cutoffs found for {protein_name} in {cutoff_file}")

    row = row.iloc[0]

    return {
        "optimized": float(row["optimized"]),
        "trusted": float(row["trusted"]),
        "noise": float(row["noise"]),
    }


def summarize_by_neighborhood(df, cutoff):
    df = df.copy()
    df["neighborhood_label"] = df["neighborhood"].map(parse_neighborhood)
    df["confusion_label"] = df.apply(lambda r: classify_row(r, cutoff), axis=1)

    summary = (
        df.groupby(["neighborhood_label", "confusion_label"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    for label in LABELS:
        if label not in summary.columns:
            summary[label] = 0

    summary["Total"] = summary[LABELS].sum(axis=1)
    summary = summary.sort_values("Total", ascending=False)

    return summary[["neighborhood_label", *LABELS, "Total"]]


def plot_summary(summary, output_path, title, top_n=30):
    if top_n is not None:
        summary = summary.head(top_n)

    summary = summary.sort_values("Total", ascending=True)

    fig_height = max(4, 0.35 * len(summary))
    fig, ax = plt.subplots(figsize=(12, fig_height))

    y = range(len(summary))
    left = [0] * len(summary)

    for label in LABELS:
        values = summary[label].values
        ax.barh(
            y,
            values,
            left=left,
            label=label,
            color=COLORS[label],
            edgecolor="black",
            linewidth=0.3,
        )
        left = [a + b for a, b in zip(left, values)]

    ax.set_yticks(y)
    ax.set_yticklabels(summary["neighborhood_label"])
    ax.set_xlabel("Number of sequences")
    ax.set_title(title)
    ax.legend(title="Classification", frameon=False)

    max_total = summary["Total"].max() if len(summary) else 0
    for i, total in enumerate(summary["Total"]):
        ax.text(total + max_total * 0.01, i, str(total), va="center", fontsize=8)

    ax.set_xlim(0, max_total * 1.12 if max_total else 1)
    fig.tight_layout()

    fig.savefig(output_path, dpi=300)
    fig.savefig(str(output_path).replace(".png", ".pdf"))
    plt.close(fig)


def process_detailed_report(report_file, cutoff_file, output_dir, top_n):
    report_file = Path(report_file)
    protein_name = report_file.name.replace("_detailed.txt", "").replace(".txt", "")

    cutoffs = load_cutoffs(cutoff_file, protein_name)
    df = pd.read_csv(report_file, sep="\t", index_col=0)

    required = {"bitscore", "training_assignment", "neighborhood"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{report_file} is missing columns: {missing}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for cutoff_name, cutoff_value in cutoffs.items():
        summary = summarize_by_neighborhood(df, cutoff_value)

        summary_file = output_dir / f"{protein_name}_{cutoff_name}_neighborhood_confusion.tsv"
        summary.to_csv(summary_file, sep="\t", index=False)

        plot_file = output_dir / f"{protein_name}_{cutoff_name}_neighborhood_confusion.png"
        title = f"{protein_name}: genomic neighbourhood confusion at {cutoff_name} cutoff"
        plot_summary(summary, plot_file, title, top_n=top_n)

        print(f"[SAVE] {summary_file}")
        print(f"[SAVE] {plot_file}")


def find_detailed_reports(path):
    path = Path(path)
    if path.is_file():
        return [path]

    return sorted(path.rglob("*_detailed.txt"))


def main():
    parser = argparse.ArgumentParser(
        description="Plot TP/FP/FN/TN composition per genomic neighbourhood from HAMSTER detailed reports."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Detailed report file or directory containing *_detailed.txt files.",
    )
    parser.add_argument(
        "-c", "--cutoffs",
        required=True,
        help="all_cutoffs.txt file with columns: name, optimized, trusted, noise.",
    )
    parser.add_argument(
        "-o", "--output",
        default="neighborhood_plots",
        help="Output directory.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=30,
        help="Show only top N neighbourhoods by total count. Use 0 for all.",
    )

    args = parser.parse_args()
    top_n = None if args.top_n == 0 else args.top_n

    reports = find_detailed_reports(args.input)
    if not reports:
        raise FileNotFoundError(f"No *_detailed.txt files found in {args.input}")

    for report in reports:
        process_detailed_report(report, args.cutoffs, args.output, top_n)


if __name__ == "__main__":
    main()