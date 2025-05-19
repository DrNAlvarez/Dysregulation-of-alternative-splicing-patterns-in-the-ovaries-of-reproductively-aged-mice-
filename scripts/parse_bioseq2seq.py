#!/usr/bin/env python3
import argparse
import sys

def parse_file(input_path):
    fasta_entries = []
    try:
        with open(input_path, 'r') as fh:
            lines = fh.readlines()
    except IOError as e:
        sys.exit(f"Error reading input file: {e}")

    i = 0
    while i < len(lines):
        id_line = lines[i].strip()
        if id_line.startswith("ID:"):
            # get the identifier
            parts = id_line.split(": ", 1)
            if len(parts) < 2:
                i += 1
                continue
            id_value = parts[1]

            # look for the next PRED line
            i += 1
            if i < len(lines):
                pred_line = lines[i].strip()
                if pred_line.startswith("PRED:"):
                    # remove leading "PRED: "
                    pred_content = pred_line[len("PRED: "):]
                    pred_vals = pred_content.split("<PC>")

                    # extract score (assumes last chunk contains "SCORE: <value>")
                    try:
                        score_value = pred_vals[-1].split("SCORE: ")[1].strip()
                    except (IndexError, ValueError):
                        score_value = ""

                    # skip non‐coding predictions
                    if "<NC>" not in pred_vals[0]:
                        # sequence is either in the 2nd chunk or the first
                        if len(pred_vals) > 1:
                            seq = pred_vals[1].split(" SCORE:")[0]
                        else:
                            seq = pred_vals[0].split(" SCORE:")[0]

                        fasta_entries.append(f">{id_value} SCORE:{score_value}")
                        fasta_entries.append(seq.strip())
        i += 1

    return fasta_entries

def write_fasta(entries, output_path):
    try:
        with open(output_path, 'w') as out:
            out.write("\n".join(entries) + "\n")
    except IOError as e:
        sys.exit(f"Error writing output file: {e}")
    print(f"Wrote {len(entries)//2} sequences to {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Parse FLAIR collapse .pep predictions into FASTA."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input predictions file"
    )
    parser.add_argument(
        "-o", "--output_prefix",
        required=True,
        help="Output file prefix (no extension)"
    )
    args = parser.parse_args()

    fasta = parse_file(args.input)
    if not fasta:
        sys.exit("No sequences found. Check your input file and parsing logic.")

    out_fp = f"{args.output_prefix}.fasta"
    write_fasta(fasta, out_fp)

if __name__ == "__main__":
    main()

