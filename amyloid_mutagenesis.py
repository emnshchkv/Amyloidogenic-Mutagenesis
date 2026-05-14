"""
Amyloidogenic Mutagenesis Script

This script performs targeted mutagenesis to reduce amyloidogenic potential by:
1. Identifying enhancer amino acids (F, Y, W, V, L, I, Q, N, G) in specified regions
2. Replacing them with beta-breaker amino acids or sequences (R, P, WY, WM, KLVFF, LPFFD)
3. Generating combinatorial mutations and outputting results to FASTA format

"""

import argparse
import sys
import os
import logging
from itertools import combinations, product
from typing import List, Tuple, Dict, Set, Optional
import re
from pathlib import Path
from datetime import datetime


def setup_logging(
    log_file: Optional[str] = None, verbose: bool = False
) -> logging.Logger:
    """
    Setup logging configuration with file output and minimal console output

    Args:
        log_file: Path to log file (required)
        verbose: If True, shows DEBUG messages in console (rarely used)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger("AmyloidMutagenesis")

    # Remove existing handlers to prevent duplicate logs
    logger.handlers = []

    # Logger accepts DEBUG and above to pass to handlers
    logger.setLevel(logging.DEBUG)

    # Create detailed formatter for file
    detailed_formatter = logging.Formatter(
        fmt="%(asctime)s - %(levelname)s - [%(funcName)s:%(lineno)d] - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Simple formatter for console (only critical info)
    simple_formatter = logging.Formatter("%(message)s")

    # Console handler - ONLY show important messages (CRITICAL, ERROR, and a few INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    if verbose:
        console_handler.setLevel(logging.DEBUG)  # Show all if verbose
    else:
        console_handler.setLevel(logging.ERROR)  # Only show errors by default
    console_handler.setFormatter(simple_formatter)
    logger.addHandler(console_handler)

    # File handler - capture EVERYTHING
    if log_file:
        try:
            file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)  # Always log DEBUG to file
            file_handler.setFormatter(detailed_formatter)
            logger.addHandler(file_handler)
        except IOError as e:
            print(f"Could not create log file '{log_file}': {e}", file=sys.stderr)
    else:
        # If no log file specified, at least show INFO to console
        console_handler.setLevel(logging.INFO)

    return logger


class AmyloidMutagenesisError(Exception):
    """Base exception for amyloid mutagenesis errors"""

    pass


class InvalidSequenceError(AmyloidMutagenesisError):
    """Exception raised for invalid protein sequences"""

    pass


class InvalidRegionError(AmyloidMutagenesisError):
    """Exception raised for invalid region specifications"""

    pass


class FASTAFileError(AmyloidMutagenesisError):
    """Exception raised for FASTA file reading errors"""

    pass


class OutputFileError(AmyloidMutagenesisError):
    """Exception raised for output file writing errors"""

    pass


class AmyloidMutagenesis:
    """Main class for performing amyloidogenic mutagenesis"""

    ENHANCERS = set("FYWVLIQNG")
    VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

    BETA_BREAKERS = {
        "single": ["R", "P"],
        "dipeptides": ["WY", "WM"],
        "pentapeptides": ["KLVFF", "LPFFD"],
    }

    def __init__(
        self,
        sequence: str,
        regions: List[Tuple[int, int]],
        output_file: str,
        logger: logging.Logger,
    ) -> None:
        """
        Initialize the mutagenesis object

        Args:
            sequence: Protein sequence string
            regions: List of tuples with (start, end) positions (1-indexed)
            output_file: Output FASTA file path
            logger: Logger instance for recording operations

        Raises:
            InvalidSequenceError: If sequence contains invalid characters
            InvalidRegionError: If regions are invalid
        """
        self.logger = logger
        self.original_sequence = sequence.upper().strip()
        self.regions = regions
        self.output_file = output_file
        self.mutations = []

        self.logger.info(
            f"Initializing AmyloidMutagenesis with sequence of length {len(self.original_sequence)}"
        )

        # Validate sequence
        self._validate_sequence()

        # Validate regions
        self._validate_regions()

        self.logger.info(
            f"Initialization successful. Sequence length: {len(self.original_sequence)}, Regions: {self.regions}"
        )

    def _validate_sequence(self) -> None:
        """
        Validate protein sequence for valid amino acids

        Raises:
            InvalidSequenceError: If sequence is empty or contains invalid characters
        """
        if not self.original_sequence:
            error_msg = "Protein sequence cannot be empty"
            self.logger.error(error_msg)
            raise InvalidSequenceError(error_msg)

        invalid_chars = set()
        for aa in self.original_sequence:
            if aa not in self.VALID_AMINO_ACIDS:
                invalid_chars.add(aa)

        if invalid_chars:
            error_msg = f"Invalid amino acid characters found: {', '.join(sorted(invalid_chars))}"
            self.logger.error(error_msg)
            raise InvalidSequenceError(error_msg)

        self.logger.debug(
            f"Sequence validation passed. Length: {len(self.original_sequence)}"
        )

    def _validate_regions(self) -> None:
        """
        Validate region specifications

        Raises:
            InvalidRegionError: If regions are invalid
        """
        if not self.regions:
            error_msg = "No regions specified"
            self.logger.error(error_msg)
            raise InvalidRegionError(error_msg)

        seq_length = len(self.original_sequence)

        for start, end in self.regions:
            if not isinstance(start, int) or not isinstance(end, int):
                error_msg = f"Region boundaries must be integers, got: {type(start).__name__}, {type(end).__name__}"
                self.logger.error(error_msg)
                raise InvalidRegionError(error_msg)

            if start < 1:
                error_msg = f"Region start position must be >= 1, got: {start}"
                self.logger.error(error_msg)
                raise InvalidRegionError(error_msg)

            if end > seq_length:
                error_msg = (
                    f"Region end position {end} exceeds sequence length {seq_length}"
                )
                self.logger.error(error_msg)
                raise InvalidRegionError(error_msg)

            if start > end:
                error_msg = f"Invalid region {start}:{end} - start position must be <= end position"
                self.logger.error(error_msg)
                raise InvalidRegionError(error_msg)

        self.logger.debug(f"Region validation passed for {len(self.regions)} regions")

    def find_enhancers_in_regions(self) -> Dict[int, str]:
        """
        Find enhancer amino acids in specified regions

        Returns:
            Dictionary mapping position (0-indexed) to amino acid
        """
        enhancers = {}

        for start, end in self.regions:
            # Convert to 0-indexed
            start_idx = start - 1
            end_idx = end - 1

            for i in range(start_idx, end_idx + 1):
                if self.original_sequence[i] in self.ENHANCERS:
                    enhancers[i] = self.original_sequence[i]
                    self.logger.debug(
                        f"Found enhancer {self.original_sequence[i]} at position {i+1}"
                    )

        if not enhancers:
            self.logger.warning("No enhancer amino acids found in specified regions")
        else:
            self.logger.info(f"Found {len(enhancers)} enhancer positions in regions")

        return enhancers

    def generate_single_mutations(
        self, enhancers: Dict[int, str]
    ) -> List[Tuple[str, str]]:
        """
        Generate single amino acid mutations

        Args:
            enhancers: Dictionary of enhancer positions and amino acids

        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        error_count = 0

        try:
            for pos, orig_aa in enhancers.items():
                for breaker in self.BETA_BREAKERS["single"]:
                    try:
                        mutant_seq = list(self.original_sequence)
                        mutant_seq[pos] = breaker
                        mutant_sequence = "".join(mutant_seq)

                        description = f"{orig_aa}{pos+1}{breaker}"
                        mutations.append((mutant_sequence, description))
                        self.logger.debug(f"Generated single mutation: {description}")
                    except Exception as e:
                        error_count += 1
                        self.logger.error(
                            f"Error generating single mutation at position {pos+1}: {e}"
                        )

            if error_count > 0:
                self.logger.warning(
                    f"Generated {len(mutations)} single mutations with {error_count} errors"
                )
            else:
                self.logger.info(
                    f"Successfully generated {len(mutations)} single mutations"
                )

        except Exception as e:
            self.logger.error(f"Fatal error in generate_single_mutations: {e}")
            raise

        return mutations

    def generate_dipeptide_replacements(
        self, enhancers: Dict[int, str]
    ) -> List[Tuple[str, str]]:
        """
        Generate dipeptide insertions at enhancer positions

        Args:
            enhancers: Dictionary of enhancer positions and amino acids

        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        error_count = 0

        try:
            for pos, orig_aa in enhancers.items():
                for dipeptide in self.BETA_BREAKERS["dipeptides"]:
                    try:
                        mutant_seq = (
                            self.original_sequence[:pos]
                            + dipeptide
                            + self.original_sequence[pos + 1 :]
                        )

                        description = f"{orig_aa}{pos+1}_{dipeptide}"
                        mutations.append((mutant_seq, description))
                        self.logger.debug(
                            f"Generated dipeptide insertion: {description}"
                        )
                    except Exception as e:
                        error_count += 1
                        self.logger.error(
                            f"Error generating dipeptide insertion at position {pos+1}: {e}"
                        )

            if error_count > 0:
                self.logger.warning(
                    f"Generated {len(mutations)} dipeptide insertions with {error_count} errors"
                )
            else:
                self.logger.info(
                    f"Successfully generated {len(mutations)} dipeptide insertions"
                )

        except Exception as e:
            self.logger.error(f"Fatal error in generate_dipeptide_replacements: {e}")
            raise

        return mutations

    def generate_pentapeptide_replacements(
        self, enhancers: Dict[int, str]
    ) -> List[Tuple[str, str]]:
        """
        Generate pentapeptide insertions at enhancer positions

        Args:
            enhancers: Dictionary of enhancer positions and amino acids

        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        error_count = 0

        try:
            for pos, orig_aa in enhancers.items():
                for pentapeptide in self.BETA_BREAKERS["pentapeptides"]:
                    try:
                        mutant_seq = (
                            self.original_sequence[:pos]
                            + pentapeptide
                            + self.original_sequence[pos + 1 :]
                        )

                        description = f"{orig_aa}{pos+1}_{pentapeptide}"
                        mutations.append((mutant_seq, description))
                        self.logger.debug(
                            f"Generated pentapeptide insertion: {description}"
                        )
                    except Exception as e:
                        error_count += 1
                        self.logger.error(
                            f"Error generating pentapeptide insertion at position {pos+1}: {e}"
                        )

            if error_count > 0:
                self.logger.warning(
                    f"Generated {len(mutations)} pentapeptide insertions with {error_count} errors"
                )
            else:
                self.logger.info(
                    f"Successfully generated {len(mutations)} pentapeptide insertions"
                )

        except Exception as e:
            self.logger.error(f"Fatal error in generate_pentapeptide_replacements: {e}")
            raise

        return mutations

    def generate_combinatorial_mutations(
        self, enhancers: Dict[int, str], max_combinations: int = 3
    ) -> List[Tuple[str, str]]:
        """
        Generate combinatorial mutations using single amino acid replacements

        Args:
            enhancers: Dictionary of enhancer positions and amino acids
            max_combinations: Maximum number of simultaneous mutations

        Returns:
            List of (mutant_sequence, description) tuples

        Raises:
            ValueError: If max_combinations is invalid
        """
        mutations = []
        error_count = 0

        try:
            if max_combinations < 2:
                error_msg = f"max_combinations must be >= 2, got: {max_combinations}"
                self.logger.error(error_msg)
                raise ValueError(error_msg)

            if len(enhancers) < 2:
                self.logger.info(
                    "Skipping combinatorial mutations: fewer than 2 enhancers found"
                )
                return mutations

            positions = list(enhancers.keys())
            self.logger.info(
                f"Generating combinatorial mutations for {len(positions)} enhancers with max_combinations={max_combinations}"
            )

            for combo_size in range(2, min(len(positions) + 1, max_combinations + 1)):
                try:
                    for pos_combo in combinations(positions, combo_size):
                        breaker_combos = list(
                            product(self.BETA_BREAKERS["single"], repeat=combo_size)
                        )

                        for breakers in breaker_combos:
                            try:
                                mutant_seq = list(self.original_sequence)
                                descriptions = []

                                for i, (pos, breaker) in enumerate(
                                    zip(pos_combo, breakers)
                                ):
                                    orig_aa = enhancers[pos]
                                    mutant_seq[pos] = breaker
                                    descriptions.append(f"{orig_aa}{pos+1}{breaker}")

                                mutant_sequence = "".join(mutant_seq)
                                description = "_".join(descriptions)
                                mutations.append((mutant_sequence, description))
                                self.logger.debug(
                                    f"Generated combinatorial mutation: {description}"
                                )
                            except Exception as e:
                                error_count += 1
                                self.logger.error(
                                    f"Error generating combinatorial mutation: {e}"
                                )

                except Exception as e:
                    error_count += 1
                    self.logger.error(
                        f"Error processing combination size {combo_size}: {e}"
                    )

            if error_count > 0:
                self.logger.warning(
                    f"Generated {len(mutations)} combinatorial mutations with {error_count} errors"
                )
            else:
                self.logger.info(
                    f"Successfully generated {len(mutations)} combinatorial mutations"
                )
            
            # Warning if combinatorial mutations exceed threshold
            if len(mutations) > 500:
                warning_msg = (
                    f"⚠️  WARNING: Generated {len(mutations)} combinatorial mutations! "
                    f"Consider reducing --max-combinations parameter to limit output size."
                )
                self.logger.warning(warning_msg)
                print(warning_msg, file=sys.stderr)

        except Exception as e:
            self.logger.error(f"Fatal error in generate_combinatorial_mutations: {e}")
            raise

        return mutations

    def generate_fixed_combinations(
        self, enhancers: Dict[int, str]
    ) -> List[Tuple[str, str]]:
        """
        Generate mutations using fixed beta-breaker combinations

        Args:
            enhancers: Dictionary of enhancer positions and amino acids

        Returns:
            List of (mutant_sequence, description) tuples
        """
        mutations = []
        error_count = 0

        try:
            if not enhancers:
                self.logger.warning("No enhancers found for fixed combinations")
                return mutations

            # Fixed combination 1: All enhancers to Proline
            try:
                mutant_seq = list(self.original_sequence)
                descriptions = []

                for pos, orig_aa in enhancers.items():
                    mutant_seq[pos] = "P"
                    descriptions.append(f"{orig_aa}{pos+1}P")

                mutant_sequence = "".join(mutant_seq)
                description = "ALL_TO_P_" + "_".join(descriptions)
                mutations.append((mutant_sequence, description))
                self.logger.debug(f"Generated fixed combination: {description}")
            except Exception as e:
                error_count += 1
                self.logger.error(f"Error generating ALL_TO_P combination: {e}")

            # Fixed combination 2: All enhancers to Arginine
            try:
                mutant_seq = list(self.original_sequence)
                descriptions = []

                for pos, orig_aa in enhancers.items():
                    mutant_seq[pos] = "R"
                    descriptions.append(f"{orig_aa}{pos+1}R")

                mutant_sequence = "".join(mutant_seq)
                description = "ALL_TO_R_" + "_".join(descriptions)
                mutations.append((mutant_sequence, description))
                self.logger.debug(f"Generated fixed combination: {description}")
            except Exception as e:
                error_count += 1
                self.logger.error(f"Error generating ALL_TO_R combination: {e}")

            # Fixed combination 3: Alternating P and R
            if len(enhancers) > 1:
                try:
                    mutant_seq = list(self.original_sequence)
                    descriptions = []
                    positions = sorted(enhancers.keys())

                    for i, pos in enumerate(positions):
                        orig_aa = enhancers[pos]
                        breaker = "P" if i % 2 == 0 else "R"
                        mutant_seq[pos] = breaker
                        descriptions.append(f"{orig_aa}{pos+1}{breaker}")

                    mutant_sequence = "".join(mutant_seq)
                    description = "ALT_P_R_" + "_".join(descriptions)
                    mutations.append((mutant_sequence, description))
                    self.logger.debug(f"Generated fixed combination: {description}")
                except Exception as e:
                    error_count += 1
                    self.logger.error(f"Error generating ALT_P_R combination: {e}")

            if error_count > 0:
                self.logger.warning(
                    f"Generated {len(mutations)} fixed combinations with {error_count} errors"
                )
            else:
                self.logger.info(
                    f"Successfully generated {len(mutations)} fixed combination mutations"
                )

        except Exception as e:
            self.logger.error(f"Fatal error in generate_fixed_combinations: {e}")
            raise

        return mutations

    def run_mutagenesis(
        self,
        include_single: bool = True,
        include_dipeptides: bool = True,
        include_pentapeptides: bool = True,
        include_combinatorial: bool = True,
        include_fixed: bool = True,
        max_combinations: int = 3,
    ) -> None:
        """
        Run the complete mutagenesis pipeline

        Args:
            include_single: Include single amino acid mutations
            include_dipeptides: Include dipeptide insertions
            include_pentapeptides: Include pentapeptide insertions
            include_combinatorial: Include combinatorial mutations
            include_fixed: Include fixed combination mutations
            max_combinations: Maximum number of simultaneous mutations for combinatorial
        """
        self.logger.info("=" * 70)
        self.logger.info("Starting Amyloidogenic Mutagenesis Pipeline")
        self.logger.info("=" * 70)
        self.logger.info(f"Original sequence: {self.original_sequence}")
        self.logger.info(f"Original sequence length: {len(self.original_sequence)}")
        self.logger.info(f"Analyzing regions: {self.regions}")

        try:
            enhancers = self.find_enhancers_in_regions()

            if not enhancers:
                self.logger.warning(
                    "No enhancer amino acids found in specified regions!"
                )
                self.logger.info("Pipeline completed with no mutations generated")
                return

            self.logger.info(f"Found {len(enhancers)} enhancer positions")
            self.logger.debug(f"Enhancer positions: {enhancers}")

            mutation_stages = []

            if include_single:
                single_muts = self.generate_single_mutations(enhancers)
                self.mutations.extend(single_muts)
                mutation_stages.append(f"Single mutations: {len(single_muts)}")

            if include_dipeptides:
                dipeptide_muts = self.generate_dipeptide_replacements(enhancers)
                self.mutations.extend(dipeptide_muts)
                mutation_stages.append(f"Dipeptide insertions: {len(dipeptide_muts)}")

            if include_pentapeptides:
                pentapeptide_muts = self.generate_pentapeptide_replacements(enhancers)
                self.mutations.extend(pentapeptide_muts)
                mutation_stages.append(
                    f"Pentapeptide insertions: {len(pentapeptide_muts)}"
                )

            if include_combinatorial and len(enhancers) > 1:
                combo_muts = self.generate_combinatorial_mutations(
                    enhancers, max_combinations
                )
                self.mutations.extend(combo_muts)
                mutation_stages.append(f"Combinatorial mutations: {len(combo_muts)}")

            if include_fixed:
                fixed_muts = self.generate_fixed_combinations(enhancers)
                self.mutations.extend(fixed_muts)
                mutation_stages.append(f"Fixed combinations: {len(fixed_muts)}")

            self.logger.info("Mutation generation summary:")
            for stage in mutation_stages:
                self.logger.info(f"  - {stage}")

            self.logger.info(f"Total mutations generated: {len(self.mutations)}")
            
            # Warning if total mutations exceed threshold
            if len(self.mutations) > 500:
                warning_msg = (
                    f"⚠️  WARNING: Generated {len(self.mutations)} total mutations! "
                    f"This may result in a large output file. Consider adjusting mutation parameters."
                )
                self.logger.warning(warning_msg)
                print(warning_msg, file=sys.stderr)
            
            self.logger.info("=" * 70)

        except Exception as e:
            self.logger.error(
                f"Fatal error during mutagenesis pipeline: {e}", exc_info=True
            )
            raise

    def write_fasta(self) -> None:
        """
        Write mutations to FASTA file

        Raises:
            OutputFileError: If writing to file fails
        """
        try:
            self.logger.info(f"Writing results to FASTA file: {self.output_file}")

            # Ensure output directory exists
            output_dir = os.path.dirname(self.output_file)
            if output_dir and not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    self.logger.info(f"Created output directory: {output_dir}")
                except OSError as e:
                    error_msg = f"Could not create output directory '{output_dir}': {e}"
                    self.logger.error(error_msg)
                    raise OutputFileError(error_msg)

            # Check if output file path is writable
            try:
                with open(self.output_file, "w") as f:
                    f.write(">Original_Sequence\n")
                    f.write(f"{self.original_sequence}\n\n")

                    for i, (sequence, description) in enumerate(self.mutations, 1):
                        f.write(f">Mutant_{i:04d}_{description}\n")
                        f.write(f"{sequence}\n\n")

                self.logger.info(
                    f"Successfully wrote {len(self.mutations)} mutations to FASTA file"
                )
                self.logger.info(f"Output file: {os.path.abspath(self.output_file)}")

            except IOError as e:
                error_msg = f"IO error while writing FASTA file: {e}"
                self.logger.error(error_msg)
                raise OutputFileError(error_msg)

        except OutputFileError:
            raise
        except Exception as e:
            error_msg = f"Unexpected error writing FASTA file: {e}"
            self.logger.error(error_msg, exc_info=True)
            raise OutputFileError(error_msg)


def parse_regions(region_string: str, logger: logging.Logger) -> List[Tuple[int, int]]:
    """
    Parse region string into list of (start, end) tuples

    Args:
        region_string: String like "1:5,15:20,30:35"
        logger: Logger instance

    Returns:
        List of (start, end) tuples

    Raises:
        InvalidRegionError: If region format is invalid
    """
    regions = []

    try:
        if not region_string or not region_string.strip():
            error_msg = "Region string cannot be empty"
            logger.error(error_msg)
            raise InvalidRegionError(error_msg)

        for region in region_string.split(","):
            region = region.strip()

            if not region:
                continue

            if ":" not in region:
                error_msg = (
                    f"Invalid region format: '{region}'. Use 'start:end' format."
                )
                logger.error(error_msg)
                raise InvalidRegionError(error_msg)

            parts = region.split(":")
            if len(parts) != 2:
                error_msg = f"Invalid region format: '{region}'. Too many colons."
                logger.error(error_msg)
                raise InvalidRegionError(error_msg)

            try:
                start, end = map(int, parts)
                regions.append((start, end))
                logger.debug(f"Parsed region: {start}:{end}")
            except ValueError as e:
                error_msg = f"Invalid region format: '{region}'. Start and end must be integers."
                logger.error(error_msg)
                raise InvalidRegionError(error_msg)

        if not regions:
            error_msg = "No valid regions found"
            logger.error(error_msg)
            raise InvalidRegionError(error_msg)

        logger.info(f"Successfully parsed {len(regions)} regions")
        return regions

    except InvalidRegionError:
        raise
    except Exception as e:
        error_msg = f"Unexpected error parsing regions: {e}"
        logger.error(error_msg)
        raise InvalidRegionError(error_msg)


def read_fasta(fasta_file: str, logger: logging.Logger) -> str:
    """
    Read sequence from FASTA file (takes first sequence)

    Args:
        fasta_file: Path to FASTA file
        logger: Logger instance

    Returns:
        Protein sequence string

    Raises:
        FASTAFileError: If file reading fails
    """
    try:
        logger.info(f"Reading FASTA file: {fasta_file}")

        # Check if file exists
        if not os.path.exists(fasta_file):
            raise FASTAFileError(f"FASTA file not found: {fasta_file}")

        # Check if path is a file
        if not os.path.isfile(fasta_file):
            raise FASTAFileError(f"Path is not a file: {fasta_file}")

        # Check file size (warn if very large)
        file_size = os.path.getsize(fasta_file)
        if file_size > 10 * 1024 * 1024:  # 10 MB
            logger.warning(f"FASTA file is large ({file_size / 1024 / 1024:.2f} MB)")

        sequence = ""
        sequence_count = 0
        line_count = 0

        try:
            with open(fasta_file, "r", encoding="utf-8") as f:
                current_header = None

                for line in f:
                    line_count += 1
                    line = line.strip()

                    if not line:  # Skip empty lines
                        continue

                    if line.startswith(">"):
                        if (
                            sequence
                        ):  # If we already have a sequence, break (take first one)
                            sequence_count += 1
                            logger.info(
                                f"Found multiple sequences in FASTA file. Using first sequence from header: '{current_header}'"
                            )
                            break
                        current_header = line[1:].split()[0]  # Get header without '>'
                        logger.debug(f"Reading sequence with header: {current_header}")
                        continue

                    sequence += line

            if not sequence:
                raise FASTAFileError(f"No sequence found in FASTA file: {fasta_file}")

            logger.info(
                f"Successfully read FASTA file: {line_count} lines, {len(sequence)} amino acids"
            )
            logger.debug(f"First 50 characters: {sequence[:50]}")

            return sequence

        except UnicodeDecodeError as e:
            raise FASTAFileError(
                f"File encoding error. The file may not be valid UTF-8: {e}"
            )

    except FASTAFileError:
        raise
    except IOError as e:
        error_msg = f"IO error while reading FASTA file: {e}"
        logger.error(error_msg)
        raise FASTAFileError(error_msg)
    except Exception as e:
        error_msg = f"Unexpected error reading FASTA file: {e}"
        logger.error(error_msg, exc_info=True)
        raise FASTAFileError(error_msg)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Amyloidogenic Mutagenesis Script - Generate mutations to reduce amyloid formation potential",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
                    Examples:
                    # Using a sequence string
                    python amyloid_mutagenesis.py --sequence "MKVLIVLLIPLASAPTVIGVK" --region "5:10,15:20" --output mutations.fasta
                    
                    # Using a FASTA file
                    python amyloid_mutagenesis.py --fasta protein.fasta --region "1:50" --output mutations.fasta
                    
                    # Custom mutation types
                    python amyloid_mutagenesis.py --fasta protein.fasta --region "1:50" --output mutations.fasta \\
                                                    --no-dipeptides --max-combinations 2
                    """,
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--sequence", "-s", type=str, help="Protein sequence as string"
    )
    input_group.add_argument(
        "--fasta", "-f", type=str, help="Path to FASTA file containing protein sequence"
    )

    # Required arguments
    parser.add_argument(
        "--region",
        "-r",
        type=str,
        required=True,
        help='Regions of interest (format: "start:end,start:end") - 1-indexed',
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Output FASTA file path"
    )

    # Optional mutation control
    parser.add_argument(
        "--no-single", action="store_true", help="Exclude single amino acid mutations"
    )
    parser.add_argument(
        "--no-dipeptides", action="store_true", help="Exclude dipeptide insertions"
    )
    parser.add_argument(
        "--no-pentapeptides",
        action="store_true",
        help="Exclude pentapeptide insertions",
    )
    parser.add_argument(
        "--no-combinatorial",
        action="store_true",
        help="Exclude combinatorial mutations",
    )
    parser.add_argument(
        "--no-fixed", action="store_true", help="Exclude fixed combination mutations"
    )
    parser.add_argument(
        "--max-combinations",
        type=int,
        default=3,
        help="Maximum number of simultaneous mutations (default: 3)",
    )

    # Logging options
    parser.add_argument(
        "--log",
        "-l",
        type=str,
        default="amyloid_mutagenesis.log",
        help="Path to log file (default: amyloid_mutagenesis.log)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging (DEBUG level on console - rarely needed)",
    )

    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(log_file=args.log, verbose=args.verbose)

    # Show startup message
    print("Amyloidogenic Mutagenesis Script started")
    print(f"Log file: {args.log}")

    logger.info(f"Script started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Command line arguments: {vars(args)}")

    exit_code = 0

    try:
        # Parse and validate regions
        regions = parse_regions(args.region, logger)

        # Read sequence
        if args.sequence:
            sequence = args.sequence
            logger.info("Using provided sequence string")
        else:
            sequence = read_fasta(args.fasta, logger)

        # Validate max_combinations argument
        if args.max_combinations < 2:
            error_msg = f"--max-combinations must be >= 2, got: {args.max_combinations}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Initialize mutagenesis object
        mutagenesis = AmyloidMutagenesis(sequence, regions, args.output, logger)

        # Run mutagenesis pipeline
        mutagenesis.run_mutagenesis(
            include_single=not args.no_single,
            include_dipeptides=not args.no_dipeptides,
            include_pentapeptides=not args.no_pentapeptides,
            include_combinatorial=not args.no_combinatorial,
            include_fixed=not args.no_fixed,
            max_combinations=args.max_combinations,
        )

        # Write output
        mutagenesis.write_fasta()

        logger.info("=" * 70)
        logger.info("Amyloidogenic mutagenesis completed successfully!")
        logger.info(f"Total mutants generated: {len(mutagenesis.mutations)}")
        logger.info(f"Results saved to: {os.path.abspath(args.output)}")
        logger.info("=" * 70)

        # Show success message
        print("Execution completed successfully")
        print(f"Generated {len(mutagenesis.mutations)} mutations")
        print(f"Results: {os.path.abspath(args.output)}")

    except InvalidSequenceError as e:
        logger.error(f"Sequence validation error: {e}", exc_info=args.verbose)
        exit_code = 2
    except InvalidRegionError as e:
        logger.error(f"Region validation error: {e}", exc_info=args.verbose)
        exit_code = 3
    except FASTAFileError as e:
        logger.error(f"FASTA file error: {e}", exc_info=args.verbose)
        exit_code = 4
    except OutputFileError as e:
        logger.error(f"Output file error: {e}", exc_info=args.verbose)
        exit_code = 5
    except ValueError as e:
        logger.error(f"Validation error: {e}", exc_info=args.verbose)
        exit_code = 6
    except AmyloidMutagenesisError as e:
        logger.error(f"Mutagenesis error: {e}", exc_info=args.verbose)
        exit_code = 7
    except KeyboardInterrupt:
        logger.error("Script interrupted by user (Ctrl+C)")
        exit_code = 130
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        exit_code = 1

    logger.info(
        f"Script finished at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} with exit code {exit_code}"
    )

    sys.exit(exit_code)


if __name__ == "__main__":
    main()