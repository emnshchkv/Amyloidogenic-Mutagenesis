import argparse
import os
import shutil
import subprocess
import sys


def parse_arguments():
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Batch process FASTA sequences using the TANGO core algorithm."
    )
    
    # Try to find Tango in system PATH or environment variable as a default
    env_tango = os.environ.get("TANGO_PATH", "")
    system_tango = shutil.which("Tango.exe") or shutil.which("Tango")
    default_tango = env_tango or system_tango

    parser.add_argument(
        "-t", "--tango",
        default=default_tango,
        help="Path to the TANGO executable (e.g., C:\\path\\to\\Tango.exe). "
             "Can be omitted if TANGO_PATH env variable is set."
    )
    parser.add_argument(
        "-f", "--fasta",
        default="../data/prp_sequences.fasta",
        help="Path to the input multi-FASTA file (default: ../data/prp_sequences.fasta)"
    )
    parser.add_argument(
        "-o", "--output",
        default="../data/prp_tango",
        help="Directory to save TANGO outputs (default: ../data/prp_tango)"
    )
    
    return parser.parse_args()


def execute_tango(tango_path, output_dir, seq_id, sequence):
    """Execute a single TANGO calculation via subprocess."""
    # Standard physiological parameters used for the amyloidogenic assay
    cmd = [
        tango_path, 
        seq_id, 
        "nt=N", 
        "ct=N", 
        "ph=7.4", 
        "te=298.15", 
        "io=0.1", 
        "tf=0", 
        f"seq={sequence}"
    ]
    
    print(f"Processing sequence: {seq_id}")
    
    try:
        # Run TANGO inside the output directory so its results are isolated there
        subprocess.run(cmd, cwd=output_dir, check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        print(f"Error executing TANGO for {seq_id}: {err.stderr.decode().strip()}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Executable not found at '{tango_path}'", file=sys.stderr)
        sys.exit(1)


def process_fasta(tango_path, fasta_path, output_dir):
    """Parse multi-FASTA file and trigger batch processing."""
    if not os.path.exists(fasta_path):
        print(f"Error: Input FASTA file not found at '{fasta_path}'", file=sys.stderr)
        sys.exit(1)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with open(fasta_path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    current_id = ""
    current_seq = ""

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_id and current_seq:
                execute_tango(tango_path, output_dir, current_id, current_seq)
            current_id = line[1:]
            current_seq = ""
        else:
            current_seq += line

    # Process the final sequence remaining in the buffer
    if current_id and current_seq:
        execute_tango(tango_path, output_dir, current_id, current_seq)


def main():
    args = parse_arguments()

    if not args.tango:
        print(
            "Error: Path to TANGO executable is required.\n"
            "Provide it via the '-t/--tango' argument or set the 'TANGO_PATH' environment variable.",
            file=sys.stderr
        )
        sys.exit(1)

    print("--- Starting TANGO Batch Prediction Pipeline ---")
    print(f"TANGO Executable: {args.tango}")
    print(f"Input FASTA:      {args.fasta}")
    print(f"Output Directory: {args.output}\n")

    process_fasta(args.tango, args.fasta, args.output)
    print("\n--- Pipeline Execution Completed Successfully ---")


if __name__ == "__main__":
    main()