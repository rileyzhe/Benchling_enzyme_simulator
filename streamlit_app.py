import csv
from pathlib import Path

import pandas as pd
import streamlit as st

from find_restriction_enzymes import (
    filter_enzymes_by_fragment_size,
    find_enzymes_by_cut_count,
    generate_gel_visualization,
    generate_plasmid_map_plotly,
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


def parse_subseq_input(text: str) -> list[tuple[str, str]]:
    """
    Parse subsequence input lines. Each line may be 'Label: SEQUENCE' or bare 'SEQUENCE'.
    Returns a list of (label, sequence) pairs with sequences uppercased and whitespace stripped.
    """
    pairs = []
    unnamed = 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if ":" in line:
            label, _, seq = line.partition(":")
            label = label.strip()
            seq = "".join(seq.split()).upper()
        else:
            unnamed += 1
            label = f"Target {unnamed}"
            seq = "".join(line.split()).upper()
        if seq:
            pairs.append((label, seq))
    return pairs


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







def render_sequence_stats(sequence: str) -> None:
    if not sequence:
        return
    st.sidebar.markdown("### Sequence summary")
    st.sidebar.write(f"**Length:** {len(sequence)} bp")


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

    sequence_input = st.text_area("Paste plasmid sequence", height=250, key="sequence_input")

    st.write("**Target subsequences** (optional)")
    _MAX_SUBSEQ = 10
    _last_filled = -1
    for _i in range(_MAX_SUBSEQ):
        if st.session_state.get(f"subseq_{_i}", "").strip():
            _last_filled = _i
    _num_boxes = min(_last_filled + 2, _MAX_SUBSEQ)
    for _i in range(_num_boxes):
        st.text_input(
            f"Target {_i + 1}",
            key=f"subseq_{_i}",
            placeholder="ATGCGT... or Label: ATGCGT...",
        )
    subseq_input = "\n".join(
        st.session_state.get(f"subseq_{_i}", "").strip()
        for _i in range(_MAX_SUBSEQ)
        if st.session_state.get(f"subseq_{_i}", "").strip()
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
        value=0,
        step=50,
    )

    max_fragment = st.number_input(
        "Maximum fragment size (bp):",
        min_value=0,
        max_value=20000,
        value=10000,
        step=50,
    )

    ladder_choice = st.selectbox(
        "Choose ladder:",
        list(LADDER_OPTIONS.keys()),
        index=0,
    )

    filter_by_subseq = st.checkbox(
        "Only show enzymes that cut inside target subsequences",
        value=False,
        help="Requires at least one subsequence above to be found in the plasmid.",
    )


    if min_fragment > max_fragment:
        st.error("Minimum fragment size must be less than or equal to maximum fragment size.")

    enzyme_subset_path = Path("joglekarlab_enzymes.txt")

    if not enzyme_subset_path.exists():
        st.error(f"{enzyme_subset_path} not found.")
        return

    with st.spinner("Loading enzyme list..."):
        enzyme_names = load_enzyme_names(str(enzyme_subset_path))

    enzyme_metadata = load_enzyme_metadata("enzyme_metadata.csv")

    st.sidebar.markdown("### Enzyme subset")
    st.sidebar.write("**Joglekar Lab Enzymes**")
    st.sidebar.write(f"{len(enzyme_names)} enzymes loaded")
    st.sidebar.write(f"{len(enzyme_metadata)} enzymes with metadata")

    # ------------------------------------------------------------------ #
    # Analyze button — runs all computation and stores results in session  #
    # state. No rendering happens here so tabs stay reactive afterwards.   #
    # ------------------------------------------------------------------ #
    if st.button("Analyze"):
        st.session_state.pop("analysis_results", None)

        if min_fragment > max_fragment:
            st.error("Minimum fragment size must be less than or equal to maximum fragment size.")
            return
        sequence = parse_sequence_input(sequence_input)
        if not sequence:
            st.error("Please paste a plasmid sequence.")
            return

        render_sequence_stats(sequence)

        # Subsequence matching
        subseq_pairs = parse_subseq_input(subseq_input)
        subseq_labels = {seq: label for label, seq in subseq_pairs}
        raw_subseqs = [seq for _, seq in subseq_pairs]
        subseq_matches = {}
        all_matched_regions = []

        if raw_subseqs:
            subseq_matches = find_subsequence_matches(sequence, raw_subseqs)
            for positions in subseq_matches.values():
                all_matched_regions.extend(positions)

        if filter_by_subseq and not all_matched_regions:
            st.error(
                "Cannot filter by subsequences: none of the entered subsequences were found in the plasmid."
                if raw_subseqs else
                "Cannot filter by subsequences: no subsequences were entered."
            )
            return

        # Enzyme search + filtering
        with st.spinner("Finding matching enzymes..."):
            enzymes_with_selected_cuts = find_enzymes_by_cut_count(
                sequence, num_cuts=num_cuts, allowed_enzymes=enzyme_names
            )
            filtered = filter_enzymes_by_fragment_size(
                enzymes_with_selected_cuts, len(sequence), min_fragment, max_fragment
            )

        if filter_by_subseq:
            filtered = [
                e for e in filtered
                if enzyme_cuts_in_regions(e[1], all_matched_regions)
            ]

        subseq_note = " and cutting inside target subsequences" if filter_by_subseq else ""

        # Subsequence display records
        subseq_display = [
            {
                "label": subseq_labels.get(subseq, subseq),
                "length": len(subseq),
                "found": bool(positions),
                "positions": positions,
            }
            for subseq, positions in subseq_matches.items()
        ]

        # Gel + table rows
        gel_bytes = None
        table_rows = []
        if filtered:
            gel_bytes = generate_gel_visualization(
                filtered, ladder_fragments=LADDER_OPTIONS[ladder_choice], output_file=None
            )
            for idx, (enzyme_name, cut_positions, fragments, _score) in enumerate(filtered, 1):
                meta = get_enzyme_meta(enzyme_metadata, enzyme_name)
                cuts_in_target = (
                    ("Yes" if enzyme_cuts_in_regions(cut_positions, all_matched_regions) else "No")
                    if all_matched_regions else "N/A"
                )
                table_rows.append({
                    "Rank": idx,
                    "Enzyme": enzyme_name,
                    "Recognition site": meta["recognition_site"],
                    "# cuts": num_cuts,
                    "Cut positions": ", ".join(str(p) for p in cut_positions),
                    "Fragments (bp)": " / ".join(str(f) for f in fragments),
                    "Optimal temp": meta["optimal_temp"],
                    "rCutSmart activity": meta["rCutSmart_activity"],
                    "Heat inactivation": meta["heat_inactivation"],
                    "Methylation sensitivity": meta["methylation_sensitivity"],
                    "Cuts in target": cuts_in_target,
                    "Notes": meta["notes"],
                })

        st.session_state["analysis_results"] = {
            "subseq_display": subseq_display,
            "enzymes_total": len(enzymes_with_selected_cuts),
            "filtered_count": len(filtered),
            "num_cuts": num_cuts,
            "min_fragment": min_fragment,
            "max_fragment": max_fragment,
            "subseq_note": subseq_note,
            "gel_bytes": gel_bytes,
            "table_rows": table_rows,
            "filter_by_subseq": filter_by_subseq,
            "sequence_length": len(sequence),
            "enzyme_cuts_all": [(e[0], e[1]) for e in filtered],
        }

    # ------------------------------------------------------------------ #
    # Results — rendered from session state so tabs and table controls     #
    # stay reactive without re-running the analysis.                       #
    # ------------------------------------------------------------------ #
    if "analysis_results" not in st.session_state:
        return

    r = st.session_state["analysis_results"]

    tabs = st.tabs(["Summary & Gel", "Enzyme Table", "Plasmid Map"])

    # ── Tab 1: Summary & Gel ──────────────────────────────────────────── #
    with tabs[0]:
        if r["subseq_display"]:
            st.subheader("Target subsequence matches")
            for item in r["subseq_display"]:
                tag = f"**{item['label']}** ({item['length']} bp)"
                if not item["found"]:
                    st.warning(f"{tag} — not found in plasmid")
                else:
                    loc = ", ".join(f"{s + 1}–{e}" for s, e in item["positions"])
                    st.success(f"{tag} — {len(item['positions'])} match(es) at: {loc}")

        st.subheader("Single enzyme results")
        st.success(f"Found {r['enzymes_total']} enzymes with exactly {r['num_cuts']} cuts.")
        st.info(
            f"{r['filtered_count']} enzymes remain after filtering fragments between "
            f"{r['min_fragment']} and {r['max_fragment']} bp{r['subseq_note']}."
        )

        if r["gel_bytes"]:
            st.image(r["gel_bytes"], use_column_width=True)
        else:
            adj = ", or adjusting the subsequence targets." if r["filter_by_subseq"] else "."
            st.warning(
                f"No enzymes passed the filter. Try lowering the minimum size, "
                f"raising the maximum size{adj}"
            )

    # ── Tab 2: Enzyme Table ───────────────────────────────────────────── #
    with tabs[1]:
        if not r["table_rows"]:
            st.info("No enzymes passed the filter.")
        else:
            col_a, col_b, col_c = st.columns([2, 2, 1])
            with col_a:
                name_filter = st.text_input(
                    "Filter by enzyme name", placeholder="e.g. EcoRI, BamHI"
                )
            with col_b:
                sort_by = st.selectbox(
                    "Sort by",
                    ["Fragment spacing (default)", "Enzyme name",
                     "rCutSmart activity", "Optimal temperature", "# cuts",
                     "First cut position"],
                )
            with col_c:
                sort_ascending = st.checkbox("Ascending", value=False)

            df = pd.DataFrame(r["table_rows"])

            if name_filter:
                df = df[df["Enzyme"].str.contains(name_filter.strip(), case=False, na=False)]

            def _parse_pct(v: str) -> float:
                try:
                    return float(v.replace("%", ""))
                except ValueError:
                    return -1.0

            def _parse_temp(v: str) -> float:
                try:
                    return float(v.replace("°C", ""))
                except ValueError:
                    return -1.0

            if sort_by == "Enzyme name":
                df = df.sort_values("Enzyme", ascending=sort_ascending)
            elif sort_by == "rCutSmart activity":
                df["_s"] = df["rCutSmart activity"].map(_parse_pct)
                df = df.sort_values("_s", ascending=sort_ascending).drop(columns=["_s"])
            elif sort_by == "Optimal temperature":
                df["_s"] = df["Optimal temp"].map(_parse_temp)
                df = df.sort_values("_s", ascending=sort_ascending).drop(columns=["_s"])
            elif sort_by == "# cuts":
                df = df.sort_values("# cuts", ascending=sort_ascending)
            elif sort_by == "First cut position":
                df["_s"] = df["Cut positions"].map(
                    lambda v: int(v.split(",")[0].strip()) if v else -1
                )
                df = df.sort_values("_s", ascending=sort_ascending).drop(columns=["_s"])

            if df.empty:
                st.info("No enzymes match the name filter.")
            else:
                def _style_rcutsmart(val: str) -> str:
                    if val == "100%":
                        return "background-color: #d4edda; color: #155724"
                    if val == "Unknown":
                        return "background-color: #e9ecef; color: #6c757d"
                    if "%" in str(val):
                        return "background-color: #fff3cd; color: #856404"
                    return ""

                def _style_heat_inact(val: str) -> str:
                    if val == "Unknown":
                        return "background-color: #e9ecef; color: #6c757d"
                    if val == "Not inactivated":
                        return "background-color: #fde8d8; color: #7d3c0a"
                    return "background-color: #d4edda; color: #155724"

                def _style_unknown(val: str) -> str:
                    return "background-color: #e9ecef; color: #6c757d" if str(val) == "Unknown" else ""

                def _style_notes(val: str) -> str:
                    return "background-color: #f8d7da; color: #721c24" if "star" in str(val).lower() else ""

                generic_cols = [
                    c for c in ("Recognition site", "Optimal temp", "Methylation sensitivity")
                    if c in df.columns
                ]
                styled = (
                    df.style
                    .map(_style_rcutsmart, subset=["rCutSmart activity"])
                    .map(_style_heat_inact, subset=["Heat inactivation"])
                    .map(_style_unknown, subset=generic_cols)
                    .map(_style_notes, subset=["Notes"])
                )
                st.dataframe(styled, use_container_width=True, hide_index=True)

    # ── Tab 3: Plasmid Map ────────────────────────────────────────────── #
    with tabs[2]:
        seq_len = r.get("sequence_length", 0)
        enzyme_cuts_all = r.get("enzyme_cuts_all", [])
        if not seq_len:
            st.info("Run analysis to generate the plasmid map.")
        else:
            enzyme_names_available = [name for name, _ in enzyme_cuts_all]
            if enzyme_names_available:
                selected_map_enzymes = st.multiselect(
                    "Select enzymes to display on map",
                    options=enzyme_names_available,
                    default=enzyme_names_available[:min(8, len(enzyme_names_available))],
                    key="plasmid_map_enzyme_select",
                )
            else:
                selected_map_enzymes = []
                st.info("No enzymes passed the filter — adjusting the fragment size range may reveal candidates.")

            enzyme_cuts_for_map = [
                (name, cuts) for name, cuts in enzyme_cuts_all if name in selected_map_enzymes
            ]
            fig = generate_plasmid_map_plotly(seq_len, r["subseq_display"], enzyme_cuts_for_map)
            st.plotly_chart(fig, use_container_width=True, config={
                "scrollZoom": False,
                "doubleClick": False,
                "displayModeBar": False,
            })

    st.sidebar.markdown("### Notes")
    st.sidebar.markdown(
        "- **Star activity** — some enzymes (e.g. EcoRI, PstI) show relaxed specificity at 50% activity in CutSmart. Use HF versions to avoid this.\n"
        "- **Cut positions** are 1-based and refer to the top strand.\n"
        "- **Fragment sizes** assume a circular plasmid."
    )


if __name__ == "__main__":
    main()
