"""
Restriction enzyme analyzer: Find enzymes with specific cut patterns and visualize virtual gel.

Usage:
    python find_restriction_enzymes.py --fasta-file plasmid.fasta [--min-fragment 500] [--max-fragment 3000]

This script analyzes a plasmid FASTA sequence for restriction enzymes that cut
exactly 3 times, evaluates them by band spacing and visibility, and generates
a static gel electrophoresis visualization. Always uses enzymes from enzyme_subset.txt.
"""

from __future__ import annotations

import sys
import argparse
from pathlib import Path
from typing import List, Tuple, Optional

from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, viridis
import base64

DEFAULT_ENZYME_LIST_FILE = "enzyme_subset.txt"


def load_sequence_from_fasta(file_path: str) -> str:
    """Load a DNA sequence from a local FASTA file."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    with path.open("r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    if not lines:
        raise ValueError(f"FASTA file is empty: {file_path}")

    if lines[0].startswith(">"):
        sequence = "".join(line for line in lines[1:] if not line.startswith(">"))
    else:
        sequence = "".join(lines)

    # Normalize: remove whitespace and convert to uppercase
    return "".join(sequence.split()).upper()


def load_enzyme_names(file_path: str) -> List[str]:
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
    allowed_enzymes: Optional[List[str]] = None,
) -> List[Tuple[str, List[int]]]:
    """
    Find restriction enzymes that cut a sequence exactly num_cuts times.

    Args:
        sequence: DNA sequence (string of ATCG characters)
        num_cuts: Number of desired cut sites
        allowed_enzymes: Optional list of enzyme names to restrict the search

    Returns:
        List of tuples: (enzyme_name, [cut_positions])
    """
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
            # Skip enzymes that fail for any reason
            continue

    return sorted(results, key=lambda x: x[0])


def calculate_fragment_sizes(
    sequence_length: int, cut_positions: List[int]
) -> List[int]:
    """
    Calculate fragment sizes from cut positions.

    For circular plasmids, fragments wrap around.

    Args:
        sequence_length: Total sequence length
        cut_positions: Sorted list of cut positions

    Returns:
        List of fragment sizes in bp (sorted descending)
    """
    if not cut_positions:
        return [sequence_length]

    # Circular plasmid: fragments between cuts and from last cut to first
    positions = sorted(cut_positions)
    fragments = []

    for i in range(len(positions)):
        if i < len(positions) - 1:
            fragment_size = positions[i + 1] - positions[i]
        else:
            # Wrap-around from last cut to first cut
            fragment_size = (sequence_length - positions[i]) + positions[0]

        if fragment_size > 0:
            fragments.append(fragment_size)

    return sorted(fragments, reverse=True)


def filter_enzymes_by_fragment_size(
    enzymes: List[Tuple[str, List[int]]],
    sequence_length: int,
    min_fragment: int,
    max_fragment: int,
) -> List[Tuple[str, List[int], List[int], float]]:
    """
    Score and rank enzymes based on fragment spacing and visibility.

    Evaluates all enzymes with 3 cuts and scores them holistically on:
    - Band spacing (prefer well-separated fragments)
    - Visibility (penalize very small bands)
    
    Args:
        enzymes: List of (enzyme_name, cut_positions) tuples
        sequence_length: Total sequence length
        min_fragment: Minimum acceptable fragment size (bp)
        max_fragment: Maximum acceptable fragment size (bp)

    Returns:
        List of tuples: (enzyme_name, cut_positions, fragment_sizes, score)
        Sorted by score (highest first)
    """
    import statistics
    
    scored = []

    for enzyme_name, cut_positions in enzymes:
        fragments = calculate_fragment_sizes(sequence_length, cut_positions)

        # Filter out bands outside acceptable range
        if not all(min_fragment <= f <= max_fragment for f in fragments):
            continue

        # Score based on spacing and visibility
        min_frag = min(fragments)
        max_frag = max(fragments)
        mean_frag = statistics.mean(fragments)
        
        # Spacing score: coefficient of variation (how spread out are the bands?)
        # Higher = more spread out = better
        stdev_frag = statistics.stdev(fragments) if len(fragments) > 1 else 0
        spacing_score = (stdev_frag / mean_frag) if mean_frag > 0 else 0
        
        # Visibility score: penalize small bands, reward adequate size
        # Bands should be at least 30% of the largest
        min_to_max_ratio = min_frag / max_frag if max_frag > 0 else 0
        visibility_score = min(min_to_max_ratio, 0.7)  # Cap at 0.7
        
        # Combined score: favor spacing but also consider visibility
        score = (spacing_score * 0.6) + (visibility_score * 0.4)
        
        scored.append((enzyme_name, cut_positions, fragments, score))

    # Sort by score (highest first)
    scored.sort(key=lambda x: x[3], reverse=True)
    return scored


def generate_gel_visualization(
    results: List[Tuple[str, List[int], List[int], float]], output_file: str = "gel.png"
) -> str:
    """
    Generate a static gel electrophoresis visualization as a PNG image.

    Fragments are displayed as horizontal bands with larger fragments near the top
    and smaller fragments migrating further down (as in real gel electrophoresis).

    Args:
        results: List of (enzyme_name, cut_positions, fragment_sizes, score)
        output_file: Path to save the image file (will be converted to .png)

    Returns:
        Path to the output PNG file
    """
    if not results:
        raise ValueError("No enzymes to visualize")

    # Convert output file to .png
    output_file = output_file.replace('.html', '.png')

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    # Get max fragment size for normalization
    all_fragments = [f for _, _, frags, _ in results for f in frags]
    max_fragment = max(all_fragments)

    # Set up lanes (one per enzyme)
    lanes = list(range(len(results)))
    lane_width = 0.7

    # Color map for bands
    cmap = viridis
    norm = Normalize(vmin=0, vmax=max_fragment)
    sm = ScalarMappable(cmap=cmap, norm=norm)

    # Draw bands for each enzyme
    for lane_idx, (enzyme_name, cut_pos, fragments, score) in enumerate(results):
        x_pos = lanes[lane_idx]

        # Sort fragments by size (descending) for better visualization
        sorted_frags = sorted(fragments, reverse=True)

        # Draw each band with spacing
        for frag_idx, fragment_size in enumerate(sorted_frags):
            # Calculate vertical position (normalized by max fragment)
            # Larger fragments at top (higher y), smaller at bottom (lower y)
            y_position = (fragment_size / max_fragment) * 90 + 5

            # Band height (proportional to log of size for better visibility)
            import math
            band_height = max(1.5, math.log(fragment_size + 1) / 2)

            # Create band rectangle
            color = sm.to_rgba(fragment_size)
            rect = mpatches.Rectangle(
                (x_pos - lane_width / 2, y_position - band_height / 2),
                lane_width,
                band_height,
                facecolor=color,
                edgecolor='black',
                linewidth=1.5,
                alpha=0.8
            )
            ax.add_patch(rect)

            # Add fragment size label
            ax.text(
                x_pos,
                y_position,
                f'{fragment_size}',
                ha='center',
                va='center',
                fontsize=9,
                fontweight='bold',
                color='white' if fragment_size > max_fragment / 2 else 'black'
            )

    # Customize plot
    ax.set_xlim(-0.5, len(results) - 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Restriction Enzyme', fontsize=12, fontweight='bold')
    ax.set_ylabel('Migration Distance (larger → smaller fragments)', fontsize=12, fontweight='bold')
    ax.set_title(
        'Virtual Gel Electrophoresis - Restriction Enzyme Digest',
        fontsize=14,
        fontweight='bold',
        pad=20
    )

    # Set x-axis ticks and labels with score if available
    enzyme_labels = [
        f"{name}\n(score: {score:.2f})" if score else name
        for name, _, _, score in results
    ]
    ax.set_xticks(lanes)
    ax.set_xticklabels(enzyme_labels, fontsize=10, fontweight='bold')

    # Add grid for readability
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Add colorbar showing fragment size scale
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Fragment Size (bp)', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Gel visualization saved to: {output_file}")
    plt.close()

    return output_file


def print_results_table(
    results: List[Tuple[str, List[int], List[int], float]]
) -> None:
    """Print a formatted table of results."""
    if not results:
        print("No enzymes found matching your criteria.")
        return

    print("\n" + "=" * 115)
    print(f"✓ Found {len(results)} enzyme(s) ranked by band spacing and visibility:\n")
    print(f"{'Rank':<6} {'Enzyme':<15} {'Fragment Sizes (bp)':<40} {'Score':<10} {'Quality':<15}")
    print("-" * 115)

    for idx, (enzyme_name, cut_pos, fragments, score) in enumerate(results, 1):
        fragments_str = " / ".join(str(f) for f in sorted(fragments, reverse=True))
        
        # Quality descriptor
        if score >= 0.5:
            quality = "Excellent"
        elif score >= 0.4:
            quality = "Good"
        elif score >= 0.3:
            quality = "Fair"
        else:
            quality = "Poor"
        
        print(f"{idx:<6} {enzyme_name:<15} {fragments_str:<40} {score:<10.3f} {quality:<15}")

    print("=" * 115)


def generate_html_report(
    results: List[Tuple[str, List[int], List[int], float]],
    fasta_file: str,
    sequence_length: int,
    enzymes_tested: int,
    enzymes_with_3_cuts: int,
    min_fragment: int,
    max_fragment: int,
    gel_png_path: str,
    report_file: str = "analysis_report.html"
) -> str:
    """
    Generate an HTML report with analysis results and gel visualization.

    Args:
        results: List of (enzyme_name, cut_positions, fragment_sizes, score)
        fasta_file: Path to the FASTA file analyzed
        sequence_length: Total sequence length
        enzymes_tested: Number of enzymes tested
        enzymes_with_3_cuts: Number of enzymes with exactly 3 cuts
        min_fragment: Minimum fragment size parameter
        max_fragment: Maximum fragment size parameter
        gel_png_path: Path to the gel.png image
        report_file: Output HTML file path

    Returns:
        Path to the generated HTML report
    """
    # Embed gel image as base64
    gel_image_b64 = ""
    if Path(gel_png_path).exists():
        with open(gel_png_path, "rb") as img_file:
            gel_image_b64 = base64.b64encode(img_file.read()).decode()

    # Build results table HTML
    results_rows = ""
    for idx, (enzyme_name, cut_pos, fragments, score) in enumerate(results, 1):
        fragments_str = " / ".join(str(f) for f in sorted(fragments, reverse=True))
        
        if score >= 0.5:
            quality = "Excellent"
        elif score >= 0.4:
            quality = "Good"
        elif score >= 0.3:
            quality = "Fair"
        else:
            quality = "Poor"
        
        results_rows += f"""
        <tr>
            <td>{idx}</td>
            <td><strong>{enzyme_name}</strong></td>
            <td>{fragments_str}</td>
            <td>{score:.3f}</td>
            <td>{quality}</td>
        </tr>
        """

    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Restriction Enzyme Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
        }}
        .summary {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .summary-item {{
            margin: 8px 0;
            font-size: 16px;
        }}
        .summary-item strong {{
            color: #2c3e50;
            min-width: 200px;
            display: inline-block;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th {{
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #bdc3c7;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        .quality-excellent {{
            background-color: #d4edda;
            color: #155724;
            padding: 4px 8px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .quality-good {{
            background-color: #d1ecf1;
            color: #0c5460;
            padding: 4px 8px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .quality-fair {{
            background-color: #fff3cd;
            color: #856404;
            padding: 4px 8px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .quality-poor {{
            background-color: #f8d7da;
            color: #721c24;
            padding: 4px 8px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .gel-container {{
            text-align: center;
            margin: 30px 0;
        }}
        .gel-container img {{
            max-width: 100%;
            height: auto;
            border: 2px solid #3498db;
            border-radius: 5px;
        }}
        .gel-title {{
            font-size: 18px;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 15px;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #bdc3c7;
            color: #7f8c8d;
            font-size: 14px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Restriction Enzyme Analysis Report</h1>
        
        <h2>Analysis Summary</h2>
        <div class="summary">
            <div class="summary-item"><strong>FASTA File:</strong> {fasta_file}</div>
            <div class="summary-item"><strong>Plasmid Length:</strong> {sequence_length:,} bp</div>
            <div class="summary-item"><strong>Enzymes Tested:</strong> {enzymes_tested} (from enzyme_subset.txt)</div>
            <div class="summary-item"><strong>Enzymes with 3 Cuts:</strong> {enzymes_with_3_cuts}</div>
            <div class="summary-item"><strong>High-Quality Enzymes Found:</strong> {len(results)}</div>
            <div class="summary-item"><strong>Fragment Size Range:</strong> {min_fragment} - {max_fragment} bp</div>
        </div>

        <h2>Enzyme Ranking Results</h2>
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Enzyme</th>
                    <th>Fragment Sizes (bp)</th>
                    <th>Score</th>
                    <th>Quality</th>
                </tr>
            </thead>
            <tbody>
                {results_rows}
            </tbody>
        </table>

        <h2>Virtual Gel Electrophoresis</h2>
        <div class="gel-container">
            <div class="gel-title">Fragment Migration Pattern (Larger → Smaller)</div>
            <img src="data:image/png;base64,{gel_image_b64}" alt="Gel Electrophoresis Visualization">
        </div>

        <div class="footer">
            <p>Report generated by Restriction Enzyme Analyzer</p>
            <p>Recommendation: Use the highest-scored enzyme for optimal gel band separation.</p>
        </div>
    </div>
</body>
</html>
"""

    # Write to file
    with open(report_file, "w") as f:
        f.write(html_content)
    
    print(f"\n✓ HTML report saved to: {report_file}")
    return report_file


def main() -> int:
    """Main entrypoint."""
    parser = argparse.ArgumentParser(
        description="Analyze plasmid restriction enzymes and visualize virtual gel."
    )
    parser.add_argument(
        "--fasta-file",
        required=True,
        help="Path to a local FASTA file containing the plasmid sequence.",
    )
    parser.add_argument(
        "--min-fragment",
        type=int,
        default=500,
        help="Minimum fragment size in bp (default: 500)",
    )
    parser.add_argument(
        "--max-fragment",
        type=int,
        default=3000,
        help="Maximum fragment size in bp (default: 3000)",
    )
    parser.add_argument(
        "--output", default="gel.png", help="Output image file for gel visualization (will be PNG)"
    )
    parser.add_argument(
        "--report", default="analysis_report.html", help="Output HTML report file (default: analysis_report.html)"
    )

    args = parser.parse_args()

    try:
        print(f"Loading sequence from FASTA file: {args.fasta_file}...")
        sequence = load_sequence_from_fasta(args.fasta_file)

        enzyme_names = load_enzyme_names(DEFAULT_ENZYME_LIST_FILE)
        print(f"Testing only {len(enzyme_names)} enzymes from {DEFAULT_ENZYME_LIST_FILE}...")

        sequence_length = len(sequence)
        print(f"Sequence length: {sequence_length} bp")

        print(f"Searching for enzymes with exactly 3 cut sites...")
        enzymes_with_3_cuts = find_enzymes_by_cut_count(
            sequence, num_cuts=3, allowed_enzymes=enzyme_names
        )
        print(f"✓ Found {len(enzymes_with_3_cuts)} enzyme(s) with 3 cuts")

        print(
            f"Evaluating all enzymes by band spacing and visibility..."
        )
        filtered = filter_enzymes_by_fragment_size(
            enzymes_with_3_cuts,
            sequence_length,
            args.min_fragment,
            args.max_fragment,
        )

        print_results_table(filtered)

        all_three_cut_results = [
            (name, cuts, calculate_fragment_sizes(sequence_length, cuts), 0.0)
            for name, cuts in enzymes_with_3_cuts
        ]

        if all_three_cut_results:
            output_file = generate_gel_visualization(all_three_cut_results, args.output)
            print(f"\nOpen {output_file} to view the gel containing all enzymes with 3 cuts.")
            print("→ Larger fragments appear near the top, smaller fragments near the bottom.")
            print(f"→ Showing {len(all_three_cut_results)} enzyme(s) in the gel.")
            print(f"→ Sorted ranking is still available in the console output and report.")

            # Generate HTML report
            report_file = generate_html_report(
                filtered,
                args.fasta_file,
                sequence_length,
                len(enzyme_names),
                len(enzymes_with_3_cuts),
                args.min_fragment,
                args.max_fragment,
                output_file,
                args.report
            )
            print(f"\n📄 Open {report_file} for a complete analysis report with visualizations.")

            if not filtered:
                print(
                    f"\n⚠️ No enzymes passed the size/visibility scoring filter."
                )
                print(
                    f"  The gel still shows all {len(all_three_cut_results)} enzymes with exactly 3 cuts."
                )
        else:
            print(
                f"\n✗ No enzymes with exactly 3 cut sites were found."
            )

        return 0

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
