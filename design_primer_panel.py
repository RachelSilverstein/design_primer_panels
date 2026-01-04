# Dependencies:
#   primer3-py==2.2.0
#   pyfaidx==0.9.0.3
#   pandas

import pandas as pd
import primer3
from pyfaidx import Fasta
import sys

# ======================================================
# CONFIGURATION — EDIT THESE VALUES
# ======================================================

GENOME_FASTA = "/Users/kleinstiverlab12/Documents/hg38.fa"

FLANK_SIZE = 200

PRODUCT_SIZE_RANGE = [150, 250]

# ---- Primary (strict) Primer3 parameters ----
PRIMER3_PARAMS_STRICT = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 23,
    "PRIMER_MAX_SIZE": 26,

    "PRIMER_MIN_TM": 63.0,
    "PRIMER_OPT_TM": 65.0,
    "PRIMER_MAX_TM": 70.0,

    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,

    "PRIMER_NUM_RETURN": 1,
}

# ---- Relaxed parameters (used only if strict fails) ----
PRIMER3_PARAMS_RELAXED = {
    **PRIMER3_PARAMS_STRICT,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 65.0,
    "PRIMER_MAX_TM": 70.0,

    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,

    "PRIMER_MAX_SIZE": 30
}

# Adapter / flap sequences (5' appended)
P5_FLAP = "AATGATACGGCGACCACCGAGATCTACAC"
P7_FLAP = "CAAGCAGAAGACGGCATACGAGAT"

# ======================================================
# FUNCTIONS
# ======================================================

def fetch_flanking_sequence(genome, chrom, start, end):
    flank_start = max(0, start - FLANK_SIZE)
    flank_end = end + FLANK_SIZE
    sequence = genome[chrom][flank_start:flank_end].seq.upper()
    return sequence, flank_start


def design_primers(template_seq, target_start_index, target_length, primer3_params):
    primer3_input = {
        "SEQUENCE_TEMPLATE": template_seq,
        "SEQUENCE_TARGET": [target_start_index, target_length],
    }

    settings = primer3_params.copy()
    settings["PRIMER_PRODUCT_SIZE_RANGE"] = [PRODUCT_SIZE_RANGE]

    return primer3.bindings.designPrimers(
        primer3_input,
        settings
    )


def process_site(genome, site_name, chrom, target_start_coord, target_end_coord):

    flanking_seq, flank_start_coord = fetch_flanking_sequence(
        genome, chrom, target_start_coord, target_end_coord
    )

    target_start_index = target_start_coord - flank_start_coord
    target_length = target_end_coord - target_start_coord

    # --------------------------------------------------
    # Attempt 1: strict parameters
    # --------------------------------------------------
    primers = design_primers(
        flanking_seq,
        target_start_index,
        target_length,
        PRIMER3_PARAMS_STRICT
    )

    status = "OK"
    message = "Designed with strict parameters"

    # --------------------------------------------------
    # Attempt 2: relaxed parameters
    # --------------------------------------------------
    if primers.get("PRIMER_PAIR_NUM_RETURNED", 0) == 0:
        primers = design_primers(
            flanking_seq,
            target_start_index,
            target_length,
            PRIMER3_PARAMS_RELAXED
        )
        status = "RELAXED"
        message = "TM/GC content/length constraints relaxed"

    results = []
    n_pairs = primers.get("PRIMER_PAIR_NUM_RETURNED", 0)

    # --------------------------------------------------
    # No primers even after relaxation
    # --------------------------------------------------
    if n_pairs == 0:
        results.append({
            "site_name": site_name,
            "chrom": chrom,
            "target_start_coord": target_start_coord,
            "target_end_coord": target_end_coord,

            "primer_pair": "NA",
            "design_status": "NO_PRIMERS",
            "design_message": "No primers found (strict + relaxed)",

            "forward_primer_no_flap": "NA",
            "reverse_primer_no_flap": "NA",
            "forward_primer_with_p5": "NA",
            "reverse_primer_with_p7": "NA",

            "forward_tm": "NA",
            "reverse_tm": "NA",
            "forward_gc": "NA",
            "reverse_gc": "NA",

            "product_size": "NA",
            "amplicon_start_coord": "NA",
            "amplicon_end_coord": "NA",
            "amplicon_length": "NA",
            "amplicon_sequence": "NA",
        })
        return results

    # --------------------------------------------------
    # Primers found
    # --------------------------------------------------
    for i in range(n_pairs):
        left_start, left_len = primers[f"PRIMER_LEFT_{i}"]
        right_start, right_len = primers[f"PRIMER_RIGHT_{i}"]

        amplicon_seq = flanking_seq[
            left_start : right_start + right_len
        ]

        amplicon_start_coord = flank_start_coord + left_start
        amplicon_end_coord = flank_start_coord + right_start + right_len

        forward_primer = primers[f"PRIMER_LEFT_{i}_SEQUENCE"]
        reverse_primer = primers[f"PRIMER_RIGHT_{i}_SEQUENCE"]

        results.append({
            "site_name": site_name,
            "chrom": chrom,
            "target_start_coord": target_start_coord,
            "target_end_coord": target_end_coord,

            "primer_pair": i + 1,
            "design_status": status,
            "design_message": message,

            "forward_primer_no_flap": forward_primer,
            "reverse_primer_no_flap": reverse_primer,
            "forward_primer_with_p5": P5_FLAP + forward_primer,
            "reverse_primer_with_p7": P7_FLAP + reverse_primer,

            "forward_tm": primers[f"PRIMER_LEFT_{i}_TM"],
            "reverse_tm": primers[f"PRIMER_RIGHT_{i}_TM"],
            "forward_gc": primers[f"PRIMER_LEFT_{i}_GC_PERCENT"],
            "reverse_gc": primers[f"PRIMER_RIGHT_{i}_GC_PERCENT"],

            "product_size": primers[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "amplicon_start_coord": amplicon_start_coord,
            "amplicon_end_coord": amplicon_end_coord,
            "amplicon_length": len(amplicon_seq),
            "amplicon_sequence": amplicon_seq,
        })

    return results


# ======================================================
# MAIN
# ======================================================

def main(input_csv, output_csv):

    genome = Fasta(GENOME_FASTA)
    df = pd.read_csv(input_csv)

    required_columns = {"site_name", "chrom", "start", "end"}
    if not required_columns.issubset(df.columns):
        sys.exit("Input CSV must contain columns: site_name, chrom, start, end")

    all_results = []

    for _, row in df.iterrows():
        try:
            all_results.extend(
                process_site(
                    genome,
                    row["site_name"],
                    row["chrom"],
                    int(row["start"]),
                    int(row["end"]),
                )
            )
        except Exception as e:
            all_results.append({
                "site_name": row["site_name"],
                "chrom": row["chrom"],
                "target_start_coord": row["start"],
                "target_end_coord": row["end"],

                "primer_pair": "NA",
                "design_status": "ERROR",
                "design_message": str(e),

                "forward_primer_no_flap": "NA",
                "reverse_primer_no_flap": "NA",
                "forward_primer_with_p5": "NA",
                "reverse_primer_with_p7": "NA",

                "forward_tm": "NA",
                "reverse_tm": "NA",
                "forward_gc": "NA",
                "reverse_gc": "NA",

                "product_size": "NA",
                "amplicon_start_coord": "NA",
                "amplicon_end_coord": "NA",
                "amplicon_length": "NA",
                "amplicon_sequence": "NA",
            })

    pd.DataFrame(all_results).to_csv(output_csv, index=False)
    print(f"Primer design complete → {output_csv}")


if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.exit(
            "Usage:\n"
            "  python design_primer_panel.py input.csv output.csv"
        )

    main(sys.argv[1], sys.argv[2])
