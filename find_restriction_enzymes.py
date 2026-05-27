from __future__ import annotations

import io
import math
import statistics
from pathlib import Path

from Bio.Restriction import AllEnzymes
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np


def build_overhang_map() -> dict[str, str]:
    """Return a dict mapping enzyme name to overhang description."""
    result = {}
    for enzyme in AllEnzymes:
        ovhg = enzyme.ovhg
        if ovhg is None:
            result[enzyme.__name__] = "Unknown"
            continue
        seq = getattr(enzyme, "ovhgseq", "") or ""
        if ovhg == 0:
            result[enzyme.__name__] = "Blunt"
        elif ovhg < 0:
            result[enzyme.__name__] = f"5' {seq}" if seq else f"5' ({abs(ovhg)} nt)"
        else:
            result[enzyme.__name__] = f"3' {seq}" if seq else f"3' ({ovhg} nt)"
    return result


def load_enzyme_names(file_path: str) -> list[str]:
    """Load a list of enzyme names from a text file."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Enzyme list file not found: {file_path}")

    with path.open("r", encoding="utf-8") as handle:
        enzymes = [line.strip() for line in handle if line.strip()]

    normalized = [name.strip() for name in enzymes if name.strip()]
    if not normalized:
        raise ValueError(f"Enzyme list file is empty: {file_path}")

    return normalized


def find_enzymes_by_cut_count(
    sequence: str,
    num_cuts: int = 3,
    allowed_enzymes: list[str] | None = None,
) -> list[tuple[str, list[int]]]:
    """Find restriction enzymes that cut a sequence exactly num_cuts times."""
    seq_obj = Seq(sequence.upper())
    allowed_set = None
    if allowed_enzymes is not None:
        allowed_set = {name.upper() for name in allowed_enzymes}

    results = []
    for enzyme in AllEnzymes:
        enzyme_name = enzyme.__name__
        if allowed_set is not None and enzyme_name.upper() not in allowed_set:
            continue
        try:
            sites = enzyme.search(seq_obj)
            if len(sites) == num_cuts:
                results.append((enzyme_name, sorted(list(sites))))
        except Exception:
            continue

    return sorted(results, key=lambda x: x[0])


def calculate_fragment_sizes(
    sequence_length: int, cut_positions: list[int]
) -> list[int]:
    """Calculate fragment sizes from cut positions for a circular plasmid."""
    if not cut_positions:
        return [sequence_length]

    positions = sorted(cut_positions)
    fragments = []

    for i in range(len(positions)):
        if i < len(positions) - 1:
            fragment_size = positions[i + 1] - positions[i]
        else:
            fragment_size = (sequence_length - positions[i]) + positions[0]

        if fragment_size > 0:
            fragments.append(fragment_size)

    return sorted(fragments, reverse=True)


def filter_enzymes_by_fragment_size(
    enzymes: list[tuple[str, list[int]]],
    sequence_length: int,
    min_fragment: int,
    max_fragment: int,
) -> list[tuple[str, list[int], list[int], float]]:
    """Score and rank enzymes based on fragment spacing and visibility."""
    scored = []

    for enzyme_name, cut_positions in enzymes:
        fragments = calculate_fragment_sizes(sequence_length, cut_positions)

        if not all(min_fragment <= f <= max_fragment for f in fragments):
            continue

        min_frag = min(fragments)
        max_frag = max(fragments)
        mean_frag = statistics.mean(fragments)

        stdev_frag = statistics.stdev(fragments) if len(fragments) > 1 else 0
        spacing_score = (stdev_frag / mean_frag) if mean_frag > 0 else 0

        min_to_max_ratio = min_frag / max_frag if max_frag > 0 else 0
        visibility_score = min(min_to_max_ratio, 0.7)

        score = (spacing_score * 0.6) + (visibility_score * 0.4)
        scored.append((enzyme_name, cut_positions, fragments, score))

    scored.sort(key=lambda x: x[3], reverse=True)
    return scored


def map_fragment_size_to_y(fragment_size: int, ladder_fragments: list[int]) -> float:
    """Map a fragment size to a gel Y position using the selected ladder."""
    ladder_sorted = sorted(ladder_fragments, reverse=True)
    if len(ladder_sorted) < 2:
        return 50.0

    max_size = ladder_sorted[0]
    min_size = ladder_sorted[-1]
    fragment_size = max(min(fragment_size, max_size), min_size)
    normalized = math.log(fragment_size / min_size) / math.log(max_size / min_size)
    return normalized * 90 + 5


def generate_gel_visualization(
    results: list[tuple[str, list[int], list[int], float]],
    ladder_fragments: list[int],
    output_file: str | None = "gel.png",
) -> str | bytes:
    """Generate a gel electrophoresis PNG. Returns bytes when output_file is None."""
    if not results:
        raise ValueError("No enzymes to visualize")

    use_buffer = output_file is None

    fig, ax = plt.subplots(figsize=(12, 8))

    lane_spacing = 1.3
    lanes = [i * lane_spacing for i in range(len(results))]
    lane_width = 0.9

    for lane_idx, (enzyme_name, cut_pos, fragments, _) in enumerate(results):
        x_pos = lanes[lane_idx]
        for fragment_size in sorted(fragments, reverse=True):
            y_position = map_fragment_size_to_y(fragment_size, ladder_fragments)
            ax.hlines(
                y=y_position,
                xmin=x_pos - lane_width / 2,
                xmax=x_pos + lane_width / 2,
                colors="#4a4a4a",
                linewidth=4.0,
                alpha=0.95,
            )

    ax.set_xlim(-0.5, lanes[-1] + 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Restriction Enzyme", fontsize=12, fontweight="bold")
    ax.set_ylabel("Migration distance (ladder-based)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Virtual Gel Electrophoresis - Restriction Enzyme Digest",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    enzyme_labels = [name for name, _, _, _ in results]
    ax.set_xticks(lanes)
    ax.set_xticklabels(enzyme_labels, fontsize=10, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    ladder_ticks = sorted(set(ladder_fragments), reverse=True)
    y_tick_positions = [map_fragment_size_to_y(s, ladder_fragments) for s in ladder_ticks]
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels([str(s) for s in ladder_ticks], fontsize=9)

    ax.yaxis.grid(True, alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)

    plt.tight_layout()

    if use_buffer:
        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", dpi=150, bbox_inches="tight", facecolor="white")
        buffer.seek(0)
        plt.close()
        return buffer.getvalue()

    plt.savefig(output_file, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    return output_file


_SUBSEQ_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
]

_ENZYME_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def generate_plasmid_map_plotly(
    sequence_length: int,
    subseq_display: list[dict],
    enzyme_cuts: list[tuple[str, list[int]]],
):
    """Return an interactive Plotly Figure of the circular plasmid map."""
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError("plotly is required for interactive plasmid maps") from exc

    def p2a(pos: int) -> float:
        return math.pi / 2 - 2 * math.pi * pos / sequence_length

    R = 1.0
    traces = []

    # Backbone circle
    theta_vals = np.linspace(0, 2 * math.pi, 720)
    traces.append(go.Scatter(
        x=np.append(R * np.cos(theta_vals), R * np.cos(theta_vals[0])),
        y=np.append(R * np.sin(theta_vals), R * np.sin(theta_vals[0])),
        mode="lines",
        line=dict(color="#333333", width=3),
        hoverinfo="skip",
        showlegend=False,
    ))

    # Tick marks + annotations
    if sequence_length <= 5_000:
        tick_interval = 500
    elif sequence_length <= 15_000:
        tick_interval = 1_000
    elif sequence_length <= 30_000:
        tick_interval = 2_000
    else:
        tick_interval = 5_000

    annotations = []
    pos = 0
    while pos < sequence_length:
        a = p2a(pos)
        traces.append(go.Scatter(
            x=[0.93 * math.cos(a), 1.07 * math.cos(a)],
            y=[0.93 * math.sin(a), 1.07 * math.sin(a)],
            mode="lines",
            line=dict(color="#777777", width=1),
            hoverinfo="skip",
            showlegend=False,
        ))
        lbl = f"{pos // 1000}k" if pos >= 1_000 else str(pos)
        annotations.append(dict(
            x=1.18 * math.cos(a), y=1.18 * math.sin(a),
            text=lbl, showarrow=False,
            font=dict(size=9, color="#555555"),
        ))
        pos += tick_interval

    # Subsequence arcs
    for si, item in enumerate(subseq_display):
        if not item["found"]:
            continue
        color = _SUBSEQ_COLORS[si % len(_SUBSEQ_COLORS)]
        first = True
        for start, end in item["positions"]:
            a_start = p2a(start)
            a_end = p2a(end)
            n_pts = max(int(abs(end - start) / sequence_length * 720), 6)
            arc = np.linspace(a_end, a_start, n_pts)
            xs = np.concatenate([1.08 * np.cos(arc), 1.18 * np.cos(arc[::-1]), [1.08 * np.cos(arc[0])]])
            ys = np.concatenate([1.08 * np.sin(arc), 1.18 * np.sin(arc[::-1]), [1.08 * np.sin(arc[0])]])
            traces.append(go.Scatter(
                x=xs, y=ys,
                fill="toself",
                fillcolor=color,
                mode="lines",
                line=dict(color=color, width=0),
                opacity=0.75,
                name=item["label"],
                legendgroup=item["label"],
                showlegend=first,
                hovertemplate=(
                    f"<b>{item['label']}</b><br>"
                    f"Position: {start + 1}–{end}<br>"
                    f"Length: {end - start} bp"
                    "<extra></extra>"
                ),
            ))
            # Label annotation at arc midpoint
            mid_a = p2a((start + end) // 2)
            annotations.append(dict(
                x=1.28 * math.cos(mid_a), y=1.28 * math.sin(mid_a),
                text=f"<b>{item['label']}</b>", showarrow=False,
                font=dict(size=9, color=color),
            ))
            first = False

    # Enzyme cut lines
    for ei, (enzyme_name, cut_positions) in enumerate(enzyme_cuts):
        color = _ENZYME_COLORS[ei % len(_ENZYME_COLORS)]
        first = True
        for cp in cut_positions:
            a = p2a(cp - 1)
            traces.append(go.Scatter(
                x=[0.78 * math.cos(a), R * math.cos(a)],
                y=[0.78 * math.sin(a), R * math.sin(a)],
                mode="lines",
                line=dict(color=color, width=2.5),
                name=enzyme_name,
                legendgroup=enzyme_name,
                showlegend=first,
                hovertemplate=(
                    f"<b>{enzyme_name}</b><br>"
                    f"Cut position: {cp}"
                    "<extra></extra>"
                ),
            ))
            first = False

    # Center text annotations
    annotations += [
        dict(x=0, y=0.10, text=f"<b>{sequence_length:,} bp</b>", showarrow=False,
             font=dict(size=14, color="#222222")),
        dict(x=0, y=-0.12, text="plasmid", showarrow=False,
             font=dict(size=10, color="#888888")),
    ]

    fig = go.Figure(data=traces)
    fig.update_layout(
        annotations=annotations,
        xaxis=dict(visible=False, range=[-1.75, 1.75], scaleanchor="y", scaleratio=1, fixedrange=True),
        yaxis=dict(visible=False, range=[-1.75, 1.75], fixedrange=True),
        dragmode=False,
        margin=dict(l=10, r=10, t=10, b=10),
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=dict(
            orientation="h",
            yanchor="top", y=-0.02,
            xanchor="center", x=0.5,
            font=dict(size=11),
        ),
        hoverlabel=dict(bgcolor="white", font_size=12),
    )
    return fig


def generate_plasmid_map(
    sequence_length: int,
    subseq_display: list[dict],
    enzyme_cuts: list[tuple[str, list[int]]],
) -> bytes:
    """Generate a circular plasmid map PNG and return as bytes."""

    def p2a(pos: int) -> float:
        """0-based position → angle (radians). 0 = 12 o'clock, clockwise."""
        return math.pi / 2 - 2 * math.pi * pos / sequence_length

    R = 1.0
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.axis("off")

    # Backbone circle
    theta = np.linspace(0, 2 * math.pi, 720)
    ax.plot(R * np.cos(theta), R * np.sin(theta), color="#333333", linewidth=3, zorder=2)

    # Tick marks and labels
    if sequence_length <= 5_000:
        tick_interval = 500
    elif sequence_length <= 15_000:
        tick_interval = 1_000
    elif sequence_length <= 30_000:
        tick_interval = 2_000
    else:
        tick_interval = 5_000

    pos = 0
    while pos < sequence_length:
        a = p2a(pos)
        ax.plot(
            [0.93 * math.cos(a), 1.07 * math.cos(a)],
            [0.93 * math.sin(a), 1.07 * math.sin(a)],
            color="#555555", linewidth=1.2, zorder=3,
        )
        lbl = f"{pos // 1000}k" if pos >= 1_000 else str(pos)
        ax.text(
            1.14 * math.cos(a), 1.14 * math.sin(a), lbl,
            ha="center", va="center", fontsize=7, color="#444444",
        )
        pos += tick_interval

    # Subsequence arcs
    found_subseqs = []
    for si, item in enumerate(subseq_display):
        if not item["found"]:
            continue
        color = _SUBSEQ_COLORS[si % len(_SUBSEQ_COLORS)]
        found_subseqs.append((item["label"], color))
        for start, end in item["positions"]:
            a_start = p2a(start)
            a_end = p2a(end)
            n_pts = max(int(abs(end - start) / sequence_length * 720), 5)
            arc = np.linspace(a_end, a_start, n_pts)
            xs = np.concatenate([1.08 * np.cos(arc), 1.18 * np.cos(arc[::-1])])
            ys = np.concatenate([1.08 * np.sin(arc), 1.18 * np.sin(arc[::-1])])
            ax.fill(xs, ys, color=color, alpha=0.75, zorder=4)
            mid_a = p2a((start + end) // 2)
            ax.text(
                1.26 * math.cos(mid_a), 1.26 * math.sin(mid_a), item["label"],
                ha="center", va="center", fontsize=7, color=color, fontweight="bold",
            )

    # Enzyme cut lines
    tab10 = plt.cm.tab10.colors
    for ei, (enzyme_name, cut_positions) in enumerate(enzyme_cuts):
        color = tab10[ei % len(tab10)]
        for cp in cut_positions:
            a = p2a(cp - 1)  # 1-based → 0-based
            ax.plot(
                [0.78 * math.cos(a), R * math.cos(a)],
                [0.78 * math.sin(a), R * math.sin(a)],
                color=color, linewidth=2.5, zorder=5, solid_capstyle="round",
            )

    # Center text
    ax.text(0, 0.08, f"{sequence_length:,} bp", ha="center", va="center",
            fontsize=12, color="#222222", fontweight="bold")
    ax.text(0, -0.10, "plasmid", ha="center", va="center",
            fontsize=9, color="#666666")

    # Legend
    handles = []
    for label, color in found_subseqs:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.75, label=label))
    for ei, (enzyme_name, _) in enumerate(enzyme_cuts):
        color = tab10[ei % len(tab10)]
        handles.append(plt.Line2D([0], [0], color=color, linewidth=2.5, label=enzyme_name))
    if handles:
        ax.legend(
            handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.10),
            ncol=min(len(handles), 5), fontsize=8, frameon=True,
        )

    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.7, 1.7)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    plt.close()
    return buf.getvalue()
