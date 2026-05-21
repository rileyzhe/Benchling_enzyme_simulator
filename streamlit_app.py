import csv
from itertools import combinations
import statistics
from pathlib import Path

import streamlit as st

from find_restriction_enzymes import (
    calculate_fragment_sizes,
    filter_enzymes_by_fragment_size,
    find_enzymes_by_cut_count,
    generate_gel_visualization,
    load_enzyme_names,
)

LADDER_OPTIONS = {
    "Life 1 kb Plus": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100],
    "Life λ DNA-HindIII": [23130, 9416, 6557, 4361, 2322, 2027, 564],
    "Life 50 bp": [5000, 2000, 1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50],
    "NEB 2-Log": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50],
    "Bioline HyperLadder 1 kb Plus": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100],
    "GeneRuler 1 kb Plus": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100],
}

METADATA_DEFAULTS = {
    "recognition_site": "Unknown",
    "optimal_temp": "Unknown",
    "buffer": "Unknown",
    "rCutSmart_activity": "Unknown",
    "heat_inactivation": "Unknown",
    "methylation_sensitivity": "Unknown",
    "notes": "",
}


@st.cache_data
def load_enzyme_metadata(path: str) -> dict:
    """
    Load enzyme metadata from a CSV file.
    Returns a dict keyed by uppercase enzyme name.
    Missing or unreadable files return an empty dict.
    """
    p = Path(path)
    if not p.exists():
        return {}
    metadata = {}
    with p.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("enzyme", "").strip()
            if name:
                metadata[name.upper()] = {
                    "recognition_site": row.get("recognition_site", "").strip() or "Unknown",
                    "optimal_temp": row.get("optimal_temp", "").strip() or "Unknown",
                    "buffer": row.get("buffer", "").strip() or "Unknown",
                    "rCutSmart_activity": row.get("rCutSmart_activity", "").strip() or "Unknown",
                    "heat_inactivation": row.get("heat_inactivation", "").strip() or "Unknown",
                    "methylation_sensitivity": row.get("methylation_sensitivity", "").strip() or "Unknown",
                    "notes": row.get("notes", "").strip(),
                }
    return metadata


def get_enzyme_meta(metadata: dict, enzyme_name: str) -> dict:
    """Return metadata for an enzyme, falling back to METADATA_DEFAULTS."""
    return metadata.get(enzyme_name.upper(), METADATA_DEFAULTS)


def parse_sequence_input(sequence_input: str) -> str:
    """Normalize raw sequence or FASTA text into a plain sequence string."""
    lines = [line.strip() for line in sequence_input.splitlines() if line.strip()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        sequence = "".join(line for line in lines[1:] if not line.startswith(">"))
    else:
        sequence = "".join(lines)
    return "".join(sequence.split()).upper()


def find_subsequence_matches(plasmid, subsequences):
    """
    Find all occurrences of each subsequence in the plasmid (forward strand).
    Returns dict: subsequence -> list of (start, end) tuples, 0-based, end exclusive.
    """
    matches = {}
    for subseq in subsequences:
        subseq_upper = subseq.upper()
        positions = []
        start = 0
        while True:
            idx = plasmid.find(subseq_upper, start)
            if idx == -1:
                break
            positions.append((idx, idx + len(subseq_upper)))
            start = idx + 1
        matches[subseq] = positions
    return matches


def enzyme_cuts_in_regions(cut_positions, regions):
    """
    cut_positions: 1-based integers (BioPython convention)
    regions: list of (start, end) tuples — 0-based, end exclusive
    """
    return any(
        start <= (cut - 1) < end
        for cut in cut_positions
        for start, end in regions
    )


def get_enzyme_subseq_coverage(cut_positions, subseq_matches):
    """Return the list of subsequences that this enzyme cuts inside."""
    return [
        subseq
        for subseq, regions in subseq_matches.items()
        if any(start <= (cut - 1) < end for cut in cut_positions for start, end in regions)
    ]


def generate_enzyme_pairs(enzymes_in_target, subseq_matches, sequence_length, max_pairs=50):
    """
    Generate ranked enzyme pairs from candidates that cut in at least one target region.
    Sorted by (coverage count desc, spacing score desc), limited to max_pairs.
    """
    results = []

    for (name1, cuts1), (name2, cuts2) in combinations(enzymes_in_target, 2):
        combined_cuts = sorted(set(cuts1) | set(cuts2))
        fragments = calculate_fragment_sizes(sequence_length, combined_cuts)

        covered1 = get_enzyme_subseq_coverage(cuts1, subseq_matches)
        covered2 = get_enzyme_subseq_coverage(cuts2, subseq_matches)
        covered_union = sorted(set(covered1 + covered2))

        if len(fragments) > 1:
            mean_f = statistics.mean(fragments)
            stdev_f = statistics.stdev(fragments)
            spacing = stdev_f / mean_f if mean_f > 0 else 0
        else:
            spacing = 0

        results.append({
            "pair": f"{name1} + {name2}",
            "covered_subseqs": covered_union,
            "coverage_count": len(covered_union),
            "combined_cuts": combined_cuts,
            "fragments": fragments,
            "spacing_score": spacing,
        })

    results.sort(key=lambda x: (-x["coverage_count"], -x["spacing_score"]))
    return results[:max_pairs]


def render_sequence_stats(sequence: str) -> None:
    if not sequence:
        return
    st.sidebar.markdown("### Sequence summary")
    st.sidebar.write(f"**Length:** {len(sequence)} bp")
    st.sidebar.write(f"**ATGC fraction:** {sum(sequence.count(nt) for nt in 'ATGC') / max(1, len(sequence)):.2%}")


def main() -> None:
    st.set_page_config(
        page_title="Restriction Enzyme Analyzer",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 Restriction Enzyme Browser")
    st.write(
        "Paste your plasmid sequence into the box below, choose the number of cuts, "
        "and filter enzymes by fragment size."
    )

    sequence_input = st.text_area("Paste plasmid sequence", height=250)

    subseq_input = st.text_area(
        "Target subsequences (optional, one per line)",
        height=100,
        help=(
            "Enter one DNA subsequence per line to define regions of interest. "
            "The app will locate each subsequence in the plasmid and can filter "
            "enzymes to only those that cut inside those regions."
        ),
    )

    num_cuts = st.number_input(
        "Desired number of cuts:",
        min_value=1,
        max_value=10,
        value=3,
        step=1,
    )

    min_fragment = st.number_input(
        "Minimum fragment size (bp):",
        min_value=0,
        max_value=10000,
        value=500,
        step=50,
    )

    max_fragment = st.number_input(
        "Maximum fragment size (bp):",
        min_value=0,
        max_value=20000,
        value=3000,
        step=50,
    )

    ladder_choice = st.selectbox(
        "Choose ladder:",
        ["Select ladder"] + list(LADDER_OPTIONS.keys()),
        index=0,
    )

    filter_by_subseq = st.checkbox(
        "Only show enzymes that cut inside target subsequences",
        value=False,
        help="Requires at least one subsequence above to be found in the plasmid.",
    )

    show_pairs = st.checkbox(
        "Show enzyme combinations that cut target subsequences",
        value=False,
        help=(
            "Generates all pairs of enzymes (from the cut-count results) that cut inside "
            "a target region. Shows up to the 50 best pairs ranked by subsequence coverage."
        ),
    )

    if min_fragment > max_fragment:
        st.error("Minimum fragment size must be less than or equal to maximum fragment size.")

    enzyme_subset_path = Path("enzyme_subset.txt")
    if not enzyme_subset_path.exists():
        st.error("enzyme_subset.txt not found. Please add your enzyme list file.")
        return

    with st.spinner("Loading enzyme subset..."):
        enzyme_names = load_enzyme_names(str(enzyme_subset_path))

    enzyme_metadata = load_enzyme_metadata("enzyme_metadata.csv")

    st.sidebar.markdown("### Enzyme subset")
    st.sidebar.write(f"{len(enzyme_names)} enzymes loaded")
    st.sidebar.write(f"{len(enzyme_metadata)} enzymes with metadata")

    if st.button("Analyze"):
        if min_fragment > max_fragment:
            st.error("Minimum fragment size must be less than or equal to maximum fragment size.")
            return

        if ladder_choice == "Select ladder":
            st.error("Please select a ladder before analyzing.")
            return

        sequence = parse_sequence_input(sequence_input)
        if not sequence:
            st.error("Please paste a plasmid sequence.")
            return

        render_sequence_stats(sequence)

        # --- Subsequence matching ---
        raw_subseqs = [line.strip() for line in subseq_input.splitlines() if line.strip()]
        subseq_matches = {}
        all_matched_regions = []

        if raw_subseqs:
            subseq_matches = find_subsequence_matches(sequence, raw_subseqs)

            st.subheader("Target subsequence matches")
            for subseq, positions in subseq_matches.items():
                if not positions:
                    st.warning(f"Subsequence not found in plasmid: `{subseq}`")
                else:
                    location_strs = [f"{start + 1}–{end}" for start, end in positions]
                    st.success(
                        f"`{subseq}` — {len(positions)} match(es) at: {', '.join(location_strs)}"
                    )
                    all_matched_regions.extend(positions)

        if filter_by_subseq and not all_matched_regions:
            if raw_subseqs:
                st.error(
                    "Cannot filter by subsequences: none of the entered subsequences were found in the plasmid."
                )
            else:
                st.error("Cannot filter by subsequences: no subsequences were entered.")
            return

        # --- Single-enzyme search ---
        with st.spinner("Finding matching enzymes..."):
            enzymes_with_selected_cuts = find_enzymes_by_cut_count(
                sequence, num_cuts=num_cuts, allowed_enzymes=enzyme_names
            )
            filtered = filter_enzymes_by_fragment_size(
                enzymes_with_selected_cuts,
                len(sequence),
                min_fragment,
                max_fragment,
            )

        if filter_by_subseq:
            filtered = [
                (name, cuts, frags, score)
                for name, cuts, frags, score in filtered
                if enzyme_cuts_in_regions(cuts, all_matched_regions)
            ]

        # --- Single-enzyme results ---
        st.subheader("Single enzyme results")

        st.success(f"Found {len(enzymes_with_selected_cuts)} enzymes with exactly {num_cuts} cuts.")
        subseq_note = " and cutting inside target subsequences" if filter_by_subseq else ""
        st.info(
            f"{len(filtered)} enzymes remain after filtering fragments between "
            f"{min_fragment} and {max_fragment} bp{subseq_note}."
        )

        if filtered:
            st.write("The gel preview below shows only enzymes that passed the fragment-size filter.")

            ladder_fragments = LADDER_OPTIONS[ladder_choice]
            gel_bytes = generate_gel_visualization(
                filtered,
                ladder_fragments=ladder_fragments,
                output_file=None,
            )
            st.image(gel_bytes, caption="Filtered enzymes with selected ladder", use_column_width=True)

            table_rows = []
            for idx, (enzyme_name, cut_positions, fragments, score) in enumerate(filtered, 1):
                meta = get_enzyme_meta(enzyme_metadata, enzyme_name)

                if all_matched_regions:
                    cuts_in_target = "Yes" if enzyme_cuts_in_regions(cut_positions, all_matched_regions) else "No"
                else:
                    cuts_in_target = "N/A"

                table_rows.append(
                    {
                        "Rank": idx,
                        "Enzyme": enzyme_name,
                        "Recognition site": meta["recognition_site"],
                        "# cuts": num_cuts,
                        "Cut positions": ", ".join(str(p) for p in cut_positions),
                        "Fragments (bp)": " / ".join(str(f) for f in fragments),
                        "Optimal temp": meta["optimal_temp"],
                        "Buffer": meta["buffer"],
                        "rCutSmart activity": meta["rCutSmart_activity"],
                        "Heat inactivation": meta["heat_inactivation"],
                        "Methylation sensitivity": meta["methylation_sensitivity"],
                        "Cuts in target": cuts_in_target,
                        "Notes": meta["notes"],
                    }
                )

            st.dataframe(table_rows, use_container_width=True)
        else:
            st.warning(
                "No enzymes passed the filter. Try lowering the minimum size, raising the maximum size"
                + (", or adjusting the subsequence targets." if filter_by_subseq else ".")
            )

        # --- Enzyme pair results ---
        if show_pairs:
            st.subheader("Enzyme pair results")

            if not all_matched_regions:
                st.info(
                    "Enter and match at least one target subsequence above to generate enzyme pairs."
                )
            else:
                enzymes_in_target = [
                    (name, cuts)
                    for name, cuts in enzymes_with_selected_cuts
                    if enzyme_cuts_in_regions(cuts, all_matched_regions)
                ]

                if len(enzymes_in_target) < 2:
                    st.info(
                        f"Only {len(enzymes_in_target)} enzyme(s) cut inside the target regions — "
                        "need at least 2 to form pairs."
                    )
                else:
                    st.write(
                        f"Generating pairs from {len(enzymes_in_target)} enzymes that cut inside "
                        "target regions. Showing up to the 50 best pairs ranked by subsequence coverage."
                    )

                    with st.spinner("Generating enzyme pairs..."):
                        pairs = generate_enzyme_pairs(
                            enzymes_in_target, subseq_matches, len(sequence), max_pairs=50
                        )

                    if pairs:
                        pair_rows = []
                        for p in pairs:
                            covered_display = ", ".join(
                                (s[:30] + "...") if len(s) > 30 else s
                                for s in p["covered_subseqs"]
                            ) or "None"
                            pair_rows.append(
                                {
                                    "Enzyme pair": p["pair"],
                                    "Targets covered": covered_display,
                                    "Combined cut positions": ", ".join(
                                        str(c) for c in p["combined_cuts"]
                                    ),
                                    "Fragment sizes (bp)": " / ".join(
                                        str(f) for f in p["fragments"]
                                    ),
                                }
                            )
                        st.dataframe(pair_rows, use_container_width=True)
                    else:
                        st.info("No pairs could be generated.")

        st.sidebar.markdown("---")
        st.sidebar.markdown("#### Notes")
        st.sidebar.write(
            "The gel image shows only enzymes that meet the selected fragment size criteria.\n"
            "Enzymes that match the cut count but fail the size filter are excluded.\n"
            "Subsequence search is performed on the forward strand only.\n"
            "Enzyme metadata sourced from enzyme_metadata.csv; 'Unknown' means no entry exists.\n"
            "Enzyme pairs are drawn from all enzymes with the selected cut count, "
            "not just those that passed the fragment-size filter."
        )


if __name__ == "__main__":
    main()
