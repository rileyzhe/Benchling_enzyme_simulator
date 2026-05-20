from pathlib import Path

import streamlit as st

from find_restriction_enzymes import (
    calculate_fragment_sizes,
    filter_enzymes_by_fragment_size,
    find_enzymes_by_cut_count,
    generate_gel_visualization,
    load_enzyme_names,
)


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

    with st.expander("How to use"):
        st.write(
            "1. Paste plasmid sequence text into the input box.\n"
            "2. Set the desired number of cuts (2 or 3 is typical).\n"
            "3. Set the minimum and maximum fragment size.\n"
            "4. Click Analyze to see the filtered enzyme list and gel preview."
        )

    sequence_input = st.text_area("Paste plasmid sequence", height=250)

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

    if min_fragment > max_fragment:
        st.error("Minimum fragment size must be less than or equal to maximum fragment size.")

    enzyme_subset_path = Path("enzyme_subset.txt")
    if not enzyme_subset_path.exists():
        st.error("enzyme_subset.txt not found. Please add your enzyme list file.")
        return

    with st.spinner("Loading enzyme subset..."):
        enzyme_names = load_enzyme_names(str(enzyme_subset_path))

    st.sidebar.markdown("### Enzyme subset")
    st.sidebar.write(f"{len(enzyme_names)} enzymes loaded")

    if st.button("Analyze"):
        if min_fragment > max_fragment:
            st.error("Minimum fragment size must be less than or equal to maximum fragment size.")
            return

        sequence = parse_sequence_input(sequence_input)
        if not sequence:
            st.error("Please paste a plasmid sequence.")
            return

        render_sequence_stats(sequence)

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

        st.success(f"Found {len(enzymes_with_selected_cuts)} enzymes with exactly {num_cuts} cuts.")
        st.info(
            f"{len(filtered)} enzymes remain after filtering fragments between {min_fragment} and {max_fragment} bp."
        )

        if filtered:
            st.subheader("Filtered enzyme gel preview")
            st.write(
                "The gel preview below shows only enzymes that passed the fragment-size filter."
            )
            gel_path = generate_gel_visualization(
                [
                    (name, cut_positions, fragments, score)
                    for name, cut_positions, fragments, score in filtered
                ],
                "streamlit_gel.png",
            )
            st.image(gel_path, caption="Filtered enzymes after size selection", use_column_width=True)

            st.subheader("Filtered enzyme list")
            table_rows = []
            for idx, (enzyme_name, cut_positions, fragments, score) in enumerate(filtered, 1):
                table_rows.append(
                    {
                        "Rank": idx,
                        "Enzyme": enzyme_name,
                        "Fragments (bp)": " / ".join(str(f) for f in fragments),
                    }
                )
            st.table(table_rows)
        else:
            st.warning(
                "No enzymes passed the fragment-size filter. Try lowering the minimum size or raising the maximum size."
            )

        st.sidebar.markdown("---")
        st.sidebar.markdown("#### Notes")
        st.sidebar.write(
            "The gel image shows only enzymes that meet the selected fragment size criteria. \n"
            "Enzymes that match the cut count but fail the size filter are excluded."
        )


if __name__ == "__main__":
    main()
