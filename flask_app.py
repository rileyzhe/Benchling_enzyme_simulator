from __future__ import annotations

import base64
import csv
import re
from pathlib import Path

from flask import Flask, jsonify, render_template, request

from find_restriction_enzymes import (
    filter_enzymes_by_fragment_size,
    find_enzymes_by_cut_count,
    generate_gel_visualization,
    generate_plasmid_map_plotly,
    load_enzyme_names,
)

app = Flask(__name__)

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
    return metadata.get(enzyme_name.upper(), METADATA_DEFAULTS)


def parse_sequence_input(sequence_input: str) -> str:
    lines = [line.strip() for line in sequence_input.splitlines() if line.strip()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        sequence = "".join(line for line in lines[1:] if not line.startswith(">"))
    else:
        sequence = "".join(lines)
    return "".join(sequence.split()).upper()


def find_subsequence_matches(plasmid: str, subsequences: list) -> dict:
    matches = {}
    for subseq in subsequences:
        subseq_upper = subseq.upper()
        positions = []
        start = 0
        while True:
            idx = plasmid.find(subseq_upper, start)
            if idx == -1:
                break
            positions.append([idx, idx + len(subseq_upper)])
            start = idx + 1
        matches[subseq] = positions
    return matches


def parse_position_range(text: str, sequence_length: int) -> tuple | None:
    """Parse '2000-3000', '2kb-3kb', '2k-3k', '1.5kb-2.5kb' into 0-based (start, end) or None."""
    m = re.match(
        r'^\s*(\d+(?:\.\d+)?)\s*(k|kb)?\s*-\s*(\d+(?:\.\d+)?)\s*(k|kb)?\s*$',
        text, re.IGNORECASE,
    )
    if not m:
        return None
    start_val = float(m.group(1)) * (1000 if m.group(2) else 1)
    end_val = float(m.group(3)) * (1000 if m.group(4) else 1)
    start = max(0, int(round(start_val)) - 1)  # 1-based → 0-based
    end = min(sequence_length, int(round(end_val)))
    if start >= end:
        return None
    return (start, end)


def parse_subseq_input(text: str) -> list:
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


def enzyme_cuts_in_regions(cut_positions: list, regions: list) -> bool:
    return any(
        start <= (cut - 1) < end
        for cut in cut_positions
        for start, end in regions
    )


def enzyme_cuts_in_all_targets(cut_positions: list, regions_by_target: list) -> bool:
    """True if the enzyme has at least one cut in every target's regions."""
    return all(
        any(
            start <= (cut - 1) < end
            for cut in cut_positions
            for (start, end) in target_regions
        )
        for target_regions in regions_by_target
    )


try:
    _ENZYME_NAMES: list[str] = load_enzyme_names("joglekarlab_enzymes.txt")
except Exception:
    _ENZYME_NAMES = []
_ENZYME_METADATA: dict = load_enzyme_metadata("enzyme_metadata.csv")


@app.route("/")
def index():
    return render_template(
        "index.html",
        ladder_options=list(LADDER_OPTIONS.keys()),
        enzyme_count=len(_ENZYME_NAMES),
        metadata_count=len(_ENZYME_METADATA),
    )


@app.route("/analyze", methods=["POST"])
def analyze():
    data = request.json

    sequence_raw = data.get("sequence", "")
    cuts_mode = data.get("cuts_mode", "exactly")
    raw_cuts = data.get("num_cuts", 3)
    if cuts_mode == "any":
        num_cuts = None
        max_cuts = None
    elif cuts_mode == "atmost":
        num_cuts = None
        max_cuts = int(raw_cuts) if raw_cuts is not None else 5
    else:
        num_cuts = int(raw_cuts) if raw_cuts is not None else 3
        max_cuts = None
    min_fragment = int(data.get("min_fragment", 0))
    max_fragment = int(data.get("max_fragment", 10000))
    ladder_choice = data.get("ladder_choice", "Life 1 kb Plus")
    subseq_input = data.get("subseq_input", "")
    filter_by_subseq = bool(data.get("filter_by_subseq", False))
    filter_all_targets = bool(data.get("filter_all_targets", False))

    if min_fragment > max_fragment:
        return jsonify({"error": "Minimum fragment size must be ≤ maximum fragment size."})

    sequence = parse_sequence_input(sequence_raw)
    if not sequence:
        return jsonify({"error": "Please paste a plasmid sequence."})

    if ladder_choice not in LADDER_OPTIONS:
        return jsonify({"error": "Invalid ladder choice."})

    subseq_pairs = parse_subseq_input(subseq_input)
    subseq_display = []
    all_matched_regions = []
    regions_by_target = []
    raw_seqs_to_find = []
    seq_label_map = {}

    for label, value in subseq_pairs:
        pos_range = parse_position_range(value, len(sequence))
        if pos_range is not None:
            start, end = pos_range
            subseq_display.append({
                "label": label, "length": end - start,
                "found": True, "positions": [[start, end]],
            })
            all_matched_regions.append([start, end])
            regions_by_target.append([[start, end]])
        else:
            raw_seqs_to_find.append(value)
            seq_label_map[value] = label

    if raw_seqs_to_find:
        seq_matches = find_subsequence_matches(sequence, raw_seqs_to_find)
        for seq, positions in seq_matches.items():
            subseq_display.append({
                "label": seq_label_map.get(seq, seq),
                "length": len(seq),
                "found": bool(positions),
                "positions": positions,
            })
            all_matched_regions.extend(positions)
            if positions:
                regions_by_target.append(positions)

    need_filter = filter_by_subseq or filter_all_targets
    if need_filter and not all_matched_regions:
        msg = (
            "Cannot filter: none of the entered subsequences or regions were found."
            if subseq_pairs
            else "Cannot filter: no subsequences or regions were entered."
        )
        return jsonify({"error": msg})

    enzymes_with_selected_cuts = find_enzymes_by_cut_count(
        sequence, num_cuts=num_cuts, max_cuts=max_cuts, allowed_enzymes=_ENZYME_NAMES
    )
    filtered = filter_enzymes_by_fragment_size(
        enzymes_with_selected_cuts, len(sequence), min_fragment, max_fragment
    )

    if filter_all_targets and regions_by_target:
        filtered = [e for e in filtered if enzyme_cuts_in_all_targets(e[1], regions_by_target)]
        subseq_note = " and cutting in every target region"
    elif filter_by_subseq and all_matched_regions:
        filtered = [e for e in filtered if enzyme_cuts_in_regions(e[1], all_matched_regions)]
        subseq_note = " and cutting in at least one target region"
    else:
        subseq_note = ""

    gel_image_b64 = None
    table_rows = []
    if filtered:
        gel_bytes = generate_gel_visualization(
            filtered, ladder_fragments=LADDER_OPTIONS[ladder_choice], output_file=None
        )
        gel_image_b64 = base64.b64encode(gel_bytes).decode("utf-8")

        for idx, (enzyme_name, cut_positions, fragments, _score) in enumerate(filtered, 1):
            meta = get_enzyme_meta(_ENZYME_METADATA, enzyme_name)
            cuts_in_target = (
                ("Yes" if enzyme_cuts_in_regions(cut_positions, all_matched_regions) else "No")
                if all_matched_regions
                else "N/A"
            )
            table_rows.append({
                "Rank": idx,
                "Enzyme": enzyme_name,
                "Recognition site": meta["recognition_site"],
                "# cuts": len(cut_positions),
                "Cut positions": ", ".join(str(p) for p in cut_positions),
                "Fragments (bp)": " / ".join(str(f) for f in fragments),
                "Optimal temp": meta["optimal_temp"],
                "rCutSmart activity": meta["rCutSmart_activity"],
                "Heat inactivation": meta["heat_inactivation"],
                "Methylation sensitivity": meta["methylation_sensitivity"],
                "Cuts in target": cuts_in_target,
                "Notes": meta["notes"],
            })

    return jsonify({
        "success": True,
        "subseq_display": subseq_display,
        "enzymes_total": len(enzymes_with_selected_cuts),
        "filtered_count": len(filtered),
        "num_cuts": num_cuts,
        "max_cuts": max_cuts,
        "cuts_mode": cuts_mode,
        "min_fragment": min_fragment,
        "max_fragment": max_fragment,
        "subseq_note": subseq_note,
        "gel_image_b64": gel_image_b64,
        "table_rows": table_rows,
        "filter_by_subseq": filter_by_subseq,
        "sequence_length": len(sequence),
        "enzyme_cuts_all": [[e[0], e[1]] for e in filtered],
    })


@app.route("/plasmid_map", methods=["POST"])
def plasmid_map_endpoint():
    data = request.json
    sequence_length = int(data.get("sequence_length", 0))
    subseq_display = data.get("subseq_display", [])
    enzyme_cuts = [(item[0], item[1]) for item in data.get("enzyme_cuts", [])]

    if not sequence_length:
        return jsonify({"error": "Missing sequence length."})

    fig = generate_plasmid_map_plotly(sequence_length, subseq_display, enzyme_cuts)
    return jsonify({"plotly_json": fig.to_json()})


if __name__ == "__main__":
    app.run(debug=True, port=5000)
