#!/usr/bin/python
import os
import glob
import subprocess
import traceback
from typing import Tuple

import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from concurrent.futures import ProcessPoolExecutor, as_completed
from matplotlib.backends.backend_pdf import PdfPages


def process_initial_plotting(options):
    output_dir = os.path.join(options.result_files_directory, "Reports")

    # Plotte die performance matrix mit den MCC Werten
    # plotting_performance(options, output_dir)

    # Plotter für die presence absence matrix histogramme
    plotting_matrix_histogram(options, output_dir)

    return


def plotting_matrix_histogram(options, tsv_dir: str):
    """
    Python-only replacement for the former R-based plotting_matrix_histogram().

    For each TSV file under tsv_dir whose name starts with "neighborhood_confusion_",
    create two multi-page A3 landscape PDFs:
      - filtered: only rows with TP>0 or FN>0
      - unfiltered: all rows

    Parallelized via ProcessPoolExecutor. Skips tasks if output PDF already exists.
    """

    # Rekursiv TSVs finden
    pattern = os.path.join(tsv_dir, "**", "neighborhood_confusion_*.tsv")
    tsv_files = glob.glob(pattern, recursive=True)

    if not tsv_files:
        print(f"No files matching 'neighborhood_confusion_*.tsv' found in {tsv_dir}.")
        return

    # Optional: fixe Achse wie im R Script (bei dir in R aktuell 500 gesetzt)
    max_x = getattr(options, "plot_confusion_max_x", 250)

    # Tasks: je TSV zwei Varianten
    tasks = []
    for tsv_file in sorted(tsv_files):
        base = os.path.splitext(tsv_file)[0]

        # Name analog zur bisherigen Logik (base_TRUE/ base_FALSE)
        # TRUE  => filter=TRUE (TP>0 or FN>0)
        # FALSE => filter=FALSE (no filter)
        for flag in ["TRUE", "FALSE"]:
            output_pdf = f"{base}_{flag}.pdf"
            if os.path.exists(output_pdf):
                print(f"[SKIP] '{output_pdf}' already exists")
            else:
                tasks.append((tsv_file, flag, output_pdf, max_x))

    if not tasks:
        print("[SKIP] All outputs already exist")
        return

    # Parallelisierung: Prozesse (robuster für matplotlib als Threads)
    # Falls du explizit options.cores hast, nutze das; sonst fallback.
    max_workers = min(len(tasks), int(getattr(options, "cores", os.cpu_count() or 1)))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(_plot_one_tsv_worker, tsv_file, flag, output_pdf, max_x): (tsv_file, flag)
            for (tsv_file, flag, output_pdf, max_x) in tasks
        }

        for future in as_completed(future_to_task):
            tsv_file, flag = future_to_task[future]
            try:
                tsv_file, flag, ok, msg = future.result()
            except Exception as exc:
                print(f"[ERROR] During the execution of '{tsv_file}' (filter={flag}) an error occurred: {exc}")

def _adjust_left_margin_for_yticklabels(fig, ax, min_left_mm=10, right_mm=10, top_bottom_in=0.5, pad_in=0.15):
    """
    Ensures y tick labels are fully inside the figure by increasing left margin dynamically.
    min_left_mm: minimum left margin in mm (keeps your baseline).
    pad_in: extra padding between labels and plot area.
    """
    # Need a renderer
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Union bbox of all y tick labels
    bboxes = [lab.get_window_extent(renderer=renderer) for lab in ax.get_yticklabels() if lab.get_text()]
    if not bboxes:
        return

    x0 = min(bb.x0 for bb in bboxes)  # leftmost pixel coordinate of labels
    fig_bbox = fig.get_window_extent(renderer=renderer)
    fig_left_px = fig_bbox.x0

    # How much labels stick out to the left (in pixels)
    overflow_px = max(0.0, fig_left_px - x0)

    # Convert px -> inches via dpi
    overflow_in = overflow_px / fig.dpi
    min_left_in = min_left_mm / 25.4
    right_in = right_mm / 25.4

    # Current margins (in figure fraction)
    cur = fig.subplotpars
    # Increase left margin by overflow + pad, but keep at least min_left_in
    new_left_in = max(min_left_in, (cur.left * fig.get_size_inches()[0]) + overflow_in + pad_in)

    # Convert inches to fraction
    new_left = new_left_in / fig.get_size_inches()[0]
    new_right = 1.0 - (right_in / fig.get_size_inches()[0])
    new_bottom = top_bottom_in / fig.get_size_inches()[1]
    new_top = 1.0 - (top_bottom_in / fig.get_size_inches()[1])

    # Clamp (avoid squeezing to zero)
    new_left = min(new_left, 0.65)
    fig.subplots_adjust(left=new_left, right=new_right, bottom=new_bottom, top=new_top)




def _plot_one_tsv_worker(
    tsv_file: str,
    flag: str,
    output_png_base: str,
    max_x: float,
) -> Tuple[str, str, bool, str]:
    """
    Worker entrypoint for multiprocessing.

    flag:
      "TRUE"  -> filtered (TP>0 or FN>0)
      "FALSE" -> unfiltered (recommended)
    Output:
      Writes PNG pages: {output_png_base}_page<N>.png
    """
    try:
        matplotlib.use("Agg", force=True)

        # --- Load TSV ---
        df = pd.read_csv(tsv_file, sep="\t")
        required = ["Neighborhood", "TP", "FP", "FN", "TN", "Total"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            return (tsv_file, flag, False, f"Missing columns {missing}")

        # --- Optional filtering (TRUE mode) ---
        if flag == "TRUE":
            df = df[(df["TP"] > 0) | (df["FN"] > 0)]

        if df.empty:
            return (tsv_file, flag, True, "No rows after filtering")

        df = df.reset_index(drop=True)

        # --- Layout constants (A3 landscape) ---
        A3_WIDTH_IN = 420 / 25.4
        A3_HEIGHT_IN = 297 / 25.4
        TOP_BOTTOM_IN = 0.5
        LEFT_MM = 10
        RIGHT_MM = 10

        # --- Bar geometry (R-equivalent proportion) ---
        bar_height_frac = 0.625  # 0.25 / (0.25 + 0.15)

        conf_order = ["TP", "FP", "FN", "TN"]
        conf_colors = {
            "TP": "#1b9e77",
            "FP": "#d95f02",
            "FN": "#7570b3",
            "TN": "#e7298a",
        }
        conf_labels = {
            "TP": "True Positives",
            "FP": "False Positives",
            "FN": "False Negatives",
            "TN": "True Negatives",
        }

        # --- Pagination ---
        category_height_in = (0.25 + 0.15) / 2.54  # 0.40 cm
        draw_area_in = A3_HEIGHT_IN - (2 * TOP_BOTTOM_IN)
        n_per_page = max(1, int(draw_area_in // category_height_in))
        n_pages = (len(df) + n_per_page - 1) // n_per_page

        def _adjust_left_margin(fig, ax, min_left_mm=10, pad_in=0.15):
            """Increase left margin so y tick labels are fully inside the canvas."""
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            bboxes = [
                lab.get_window_extent(renderer=renderer)
                for lab in ax.get_yticklabels()
                if lab.get_text()
            ]
            if not bboxes:
                return

            x0 = min(bb.x0 for bb in bboxes)
            fig_bbox = fig.get_window_extent(renderer=renderer)
            overflow_px = max(0.0, fig_bbox.x0 - x0)
            overflow_in = overflow_px / fig.dpi

            min_left_in = min_left_mm / 25.4
            cur_left_in = fig.subplotpars.left * fig.get_size_inches()[0]
            new_left_in = max(min_left_in, cur_left_in + overflow_in + pad_in)

            # Clamp to keep some plotting area
            new_left = min(new_left_in / fig.get_size_inches()[0], 0.65)

            fig.subplots_adjust(
                left=new_left,
                right=1.0 - (RIGHT_MM / 25.4) / fig.get_size_inches()[0],
                bottom=TOP_BOTTOM_IN / fig.get_size_inches()[1],
                top=fig.subplotpars.top,  # keep current top (we adjust below for legend)
            )

        # --- Plot pages ---
        for page_idx in range(n_pages):
            start = page_idx * n_per_page
            end = min(len(df), (page_idx + 1) * n_per_page)
            d = df.iloc[start:end].copy()

            neighborhoods = d["Neighborhood"].tolist()
            n = len(neighborhoods)
            y = list(range(n))

            fig = plt.figure(figsize=(A3_WIDTH_IN, A3_HEIGHT_IN))
            ax = fig.add_subplot(111)

            # Baseline margins. Reserve extra space at top for legend (no overlap).
            # We'll tighten further after we place legend.
            fig.subplots_adjust(
                left=(LEFT_MM / 25.4) / A3_WIDTH_IN,
                right=1.0 - (RIGHT_MM / 25.4) / A3_WIDTH_IN,
                bottom=TOP_BOTTOM_IN / A3_HEIGHT_IN,
                top=1.0 - (TOP_BOTTOM_IN / A3_HEIGHT_IN) - 0.06,  # extra headroom for legend
            )

            # --- Stacked horizontal bars (clipped at max_x) ---
            lefts = [0.0] * n
            for conf in conf_order:
                vals = d[conf].astype(float).tolist()
                widths = []
                for i, v in enumerate(vals):
                    remaining = max_x - lefts[i]
                    widths.append(0.0 if remaining <= 0 else min(v, remaining))

                ax.barh(
                    y,
                    widths,
                    left=lefts,
                    height=bar_height_frac,
                    color=conf_colors[conf],
                    edgecolor=None,
                    label=conf_labels[conf],  # expanded label
                )
                lefts = [l + w for l, w in zip(lefts, widths)]

            # --- Total labels ---
            offset = max_x * 0.01
            for i, tot in enumerate(d["Total"].astype(float)):
                x = min(tot, max_x) + offset
                ax.text(
                    x, y[i],
                    f"{int(tot) if float(tot).is_integer() else tot}",
                    va="center",
                    ha="left",
                    fontsize=12,
                    clip_on=False,
                )

            # --- Y axis ---
            ax.set_yticks(y)
            ax.set_yticklabels(neighborhoods, fontsize=12)
            for lab in ax.get_yticklabels():
                lab.set_horizontalalignment("right")
            ax.invert_yaxis()

            # Fix left margin to prevent y-label clipping
            _adjust_left_margin(fig, ax, min_left_mm=LEFT_MM)

            # --- X axis ---
            ax.set_xlim(0, max_x * 1.10)
            ax.set_xlabel("Absolute Count", fontsize=12)

            # --- No title (explicitly) ---
            ax.set_title("")

            # --- Legend: outside plot area, no overlap ---
            # Put legend centered ABOVE the axes, in the reserved top margin.
            if page_idx == 0:
                leg = ax.legend(
                    loc="lower center",
                    bbox_to_anchor=(0.5, 1.02),  # just above axes
                    ncol=2,                      # 2 columns fits longer labels better
                    frameon=False,
                    borderaxespad=0.0,
                    fontsize=12,
                )
                # Ensure enough top margin for the legend (in case it is taller than expected)
                fig.canvas.draw()
                bb = leg.get_window_extent(fig.canvas.get_renderer())
                leg_height_in = bb.height / fig.dpi
                # add a small buffer
                needed_top = (TOP_BOTTOM_IN + leg_height_in + 0.15) / fig.get_size_inches()[1]
                # Convert to top fraction: 1 - needed_top
                fig.subplots_adjust(top=1.0 - needed_top)

            out_png = f"{output_png_base}_page{page_idx + 1}.png"
            fig.savefig(out_png, dpi=300, bbox_inches="tight")
            plt.close(fig)

        return (tsv_file, flag, True, "OK")

    except Exception as e:
        return (tsv_file, flag, False, f"{e}\n{traceback.format_exc()}")
