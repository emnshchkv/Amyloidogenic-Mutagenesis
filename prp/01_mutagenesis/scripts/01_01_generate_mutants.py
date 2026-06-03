import argparse
import os
import sys


def read_fasta(file_path):
    """
    Read a single wild-type sequence from a FASTA file.
    Ignores header lines starting with '>'.
    """
    if not os.path.exists(file_path):
        print(f"Error: Input file not found at '{file_path}'", file=sys.stderr)
        sys.exit(1)

    sequence = ""
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">") and line:
                sequence += line
    return sequence


def generate_mutants(sequence, start_res, end_res, breakers):
    """
    Generate single point mutant sequences.

    :param sequence: Wild-type protein sequence (string).
    :param start_res: Start coordinate of the target region (1-based index).
    :param end_res: End coordinate of the target region (inclusive).
    :param breakers: List of amino acids to introduce as mutations.
    :return: List of dictionaries containing mutation IDs and sequences.
    """
    mutants = []
    # Convert biological 1-based numbering to Python 0-based indexing
    start_idx = start_res - 1
    end_idx = end_res

    # Validate coordinates
    if start_idx < 0 or end_idx > len(sequence) or start_idx >= end_idx:
        print(
            f"Error: Invalid coordinates. Sequence length is {len(sequence)}.",
            file=sys.stderr,
        )
        sys.exit(1)

    for i in range(start_idx, end_idx):
        original_aa = sequence[i]
        for b in breakers:
            if b == original_aa:
                continue  # Skip if the amino acid is already the target one

            # Create mutation
            mut_seq = sequence[:i] + b + sequence[i + 1 :]
            mut_name = f"{original_aa}{i + 1}{b}"

            mutants.append({"id": mut_name, "seq": mut_seq})

    return mutants


def main():
    parser = argparse.ArgumentParser(
        description="Generate single point mutations (breakers) for a target protein region."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input FASTA file containing the wild-type sequence.",
    )
    parser.add_argument(
        "-s",
        "--start",
        type=int,
        required=True,
        help="Start coordinate of the target region (1-based index).",
    )
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        required=True,
        help="End coordinate of the target region (inclusive, 1-based index).",
    )
    parser.add_argument(
        "-b",
        "--breakers",
        nargs="+",
        default=["P", "D", "E", "K", "R"],
        help="List of amino acids to substitute (default: P D E K R). Separate with spaces.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="prp_sequences.fasta",
        help="Name of the output FASTA file (default: prp_sequences.fasta).",
    )

    args = parser.parse_args()

    print("--- Starting Mutagenesis Pipeline ---")
    print(f"Input file: {args.input}")
    print(f"Target region: {args.start}-{args.end}")
    print(f"Breaker amino acids: {args.breakers}")

    # Extract WT sequence
    wt_seq = read_fasta(args.input)
    print(f"Wild-type sequence loaded. Length: {len(wt_seq)} amino acids.")

    # Generate mutants
    mutants = generate_mutants(wt_seq, args.start, args.end, args.breakers)
    print(f"Generated {len(mutants)} mutant sequences.")

    # Determine output path (same directory as input file)
    input_dir = os.path.dirname(os.path.abspath(args.input))
    output_path = os.path.join(input_dir, args.output)

    # Write to FASTA
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            # Write wild-type for comparison baseline
            f.write(f">WT\n{wt_seq}\n")
            # Write all mutations
            for m in mutants:
                f.write(f">{m['id']}\n{m['seq']}\n")
        print(f"Successfully saved all sequences to: {output_path}")
    except IOError as err:
        print(f"Error writing to output file: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
