from __future__ import annotations

import io
import math
import statistics
from pathlib import Path

from Bio.Restriction import AllEnzymes
from Bio.Seq import Seq
import matplotlib.pyplot as plt


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
