import argparse
import os
import sys


def trim_sequences(input_fasta, output_fasta, start_res, end_res):
    """
    Trim FASTA sequences to the specified biological coordinates.

    :param input_fasta: Path to the input FASTA file.
    :param output_fasta: Path to save the trimmed FASTA file.
    :param start_res: 1-based start coordinate of the sequence to keep.
    :param end_res: 1-based end coordinate of the sequence to keep.
    :return: Number of processed sequences.
    """
    processed_count = 0
    # Convert biological 1-based start coordinate to Python 0-based index
    start_idx = start_res - 1
    # End coordinate remains the same since Python slice [a:b] is exclusive at b
    end_idx = end_res

    try:
        with open(input_fasta, "r", encoding="utf-8") as infile, \
                open(output_fasta, "w", encoding="utf-8") as outfile:

            header = ""
            sequence = []

            def write_record():
                nonlocal processed_count
                full_seq = "".join(sequence)

                # Check if sequence is long enough to be trimmed
                if len(full_seq) < end_idx:
                    print(
                        f"Warning: Sequence {header[1:]} is shorter than end coordinate ({len(full_seq)} < {end_idx}).",
                        file=sys.stderr)

                trimmed_seq = full_seq[start_idx:end_idx]

                outfile.write(f"{header}\n{trimmed_seq}\n")
                processed_count += 1

            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if header:
                        write_record()
                    header = line
                    sequence = []
                else:
                    sequence.append(line)

            # Do not forget to write the last sequence in the buffer
            if header:
                write_record()

    except IOError as err:
        print(f"Error processing files: {err}", file=sys.stderr)
        sys.exit(1)

    return processed_count


def main():
    parser = argparse.ArgumentParser(
        description="Trim N- and C-terminal regions from multi-FASTA sequences."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input FASTA/multi-FASTA file."
    )
    parser.add_argument(
        "-s", "--start",
        type=int,
        required=True,
        help="Start coordinate of the region to keep (1-based index, e.g., 23)."
    )
    parser.add_argument(
        "-e", "--end",
        type=int,
        required=True,
        help="End coordinate of the region to keep (1-based index, e.g., 231)."
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file not found at '{args.input}'", file=sys.stderr)
        sys.exit(1)

    if args.start >= args.end or args.start < 1:
        print("Error: Invalid coordinates. Ensure that start < end and start >= 1.", file=sys.stderr)
        sys.exit(1)

    # Generate output filename dynamically
    input_dir = os.path.dirname(os.path.abspath(args.input))
    input_basename = os.path.basename(args.input)
    file_name, file_ext = os.path.splitext(input_basename)

    output_filename = f"{file_name}_{args.start}-to-{args.end}{file_ext}"
    output_path = os.path.join(input_dir, output_filename)

    print("--- Starting Sequence Trimming ---")
    print(f"Input file:   {args.input}")
    print(f"Keep region:  {args.start} to {args.end}")

    # Run the trimming process
    count = trim_sequences(args.input, output_path, args.start, args.end)

    print(f"Output file:  {output_path}")
    print(f"Successfully processed {count} sequences.")


if __name__ == "__main__":
    main()