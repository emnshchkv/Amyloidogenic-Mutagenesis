"""
Single Amino Acid Substitution Generator for Protein Sequences

This script generates all possible single amino acid mutants within specified
regions of a protein sequence of Abeta39 and outputs them in FASTA format.

The script includes the wild-type sequence and all single-point mutants
(excluding the original amino acid at each position).
"""

import argparse
from typing import List, Tuple

# All 20 standard amino acids (one-letter codes)
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def find_region_positions(sequence: str, region: str) -> List[int]:
    """
    Find the start position of a region within the sequence and return
    all 0-based indices occupied by that region.

    Args:
        sequence (str): Full protein sequence.
        region (str): Subsequence (region) to locate.

    Returns:
        List[int]: List of positions (0-based) corresponding to the region.

    Raises:
        ValueError: If the region is not found in the sequence.
    """
    start = sequence.find(region)
    if start == -1:
        raise ValueError(f"Region '{region}' not found in the sequence.")

    return list(range(start, start + len(region)))


def generate_mutants(sequence: str, positions: List[int]) -> List[Tuple[str, str]]:
    """
    Generate all single amino acid mutants for the given positions.

    Args:
        sequence (str): Original protein sequence.
        positions (List[int]): List of 0-based positions to mutate.

    Returns:
        List[Tuple[str, str]]: List of tuples containing (header, mutated_sequence).
    """
    mutants = []

    for pos in positions:
        original_aa = sequence[pos]

        for new_aa in AMINO_ACIDS:
            if new_aa == original_aa:
                continue

            # Create mutated sequence
            mutated = sequence[:pos] + new_aa + sequence[pos + 1 :]

            # Create header in standard mutation notation (e.g., K21M)
            header = f">mutant_{original_aa}{pos + 1}{new_aa}"
            mutants.append((header, mutated))

    return mutants


def main() -> None:
    """Main function to parse arguments and generate mutant library."""

    parser = argparse.ArgumentParser(
        description="Generate all single amino acid substitutions within specified protein regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--seq",
        type=str,
        default="DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGV",
        help="Full amino acid sequence of the wild-type protein.",
    )

    parser.add_argument(
        "--regions",
        type=str,
        nargs="+",
        default=["QKLVFFA", "SNKGAIIGLMVGGV"],
        help="One or more regions (subsequences) to mutate. "
        "Mutations will be generated at every position within these regions.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="mutants.fasta",
        help="Output FASTA filename.",
    )

    args, unknown = parser.parse_known_args()

    if unknown:
        print(f"Warning: Ignoring unknown arguments: {unknown}")

    sequence = args.seq
    regions = args.regions
    output_file = args.output

    # Collect all positions to mutate
    positions: List[int] = []
    for region in regions:
        try:
            positions.extend(find_region_positions(sequence, region))
        except ValueError as e:
            print(f"Error: {e}")
            return

    # Remove duplicate positions (in case of overlapping regions) and sort
    positions = sorted(set(positions))

    print(f"Found {len(positions)} unique positions for mutation.")

    # Generate mutants
    mutants = generate_mutants(sequence, positions)

    # Validation
    expected_mutants = len(positions) * 19
    if len(mutants) != expected_mutants:
        print(
            f"Warning: Generated {len(mutants)} mutants, expected {expected_mutants}."
        )

    # Write to FASTA file
    with open(output_file, "w") as f:
        # Write wild-type sequence first
        f.write(">WT\n")
        f.write(sequence + "\n")

        # Write all mutants
        for header, mutated_seq in mutants:
            f.write(header + "\n")
            f.write(mutated_seq + "\n")

    total_entries = 1 + len(mutants)
    print(
        f"Mutant library successfully written to '{output_file}' "
        f"({total_entries:,} total entries: 1 WT + {len(mutants):,} mutants)."
    )


if __name__ == "__main__":
    main()
