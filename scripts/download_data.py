#!/usr/bin/env python3
"""
Download heavy data assets for the Amyloidogenic-Mutagenesis project.

Large binaries (trajectories, structures, force field, animations) are kept out
of git to keep the repository clonable and container-friendly. They live in a
shared Google Drive folder and are fetched on demand by this script.

Usage
-----
    pip install gdown
    python scripts/download_data.py                 # fetch everything
    python scripts/download_data.py --only prp       # one subproject
    python scripts/download_data.py --list           # show what would be fetched
    python scripts/download_data.py --dest /tmp/data # custom staging location

Layout expectation
-------------------
The Drive folder should mirror the repository's data layout, i.e. it should
contain top-level folders `abeta42/`, `abeta39/`, `prp/`, and `shared/`
(the force field). gdown downloads the folder tree as-is into the repository
root, so files land in the directories the code already expects.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Shared Google Drive folder holding all heavy assets.
DRIVE_FOLDER_URL = (
    "https://drive.google.com/drive/folders/"
    "1pGoWKy3sahO6ZpG4wmWFTG1DNPLJWU4M?usp=share_link"
)

# Subproject -> Drive subfolder name. Keep these in sync with the Drive layout.
SUBPROJECTS = {
    "abeta42": "abeta42",
    "abeta39": "abeta39",
    "prp": "prp",
    "shared": "shared",  # CHARMM36m force field, shared across all three
}

REPO_ROOT = Path(__file__).resolve().parent.parent


def check_gdown() -> None:
    if shutil.which("gdown") is None:
        sys.exit(
            "gdown is not installed.\n  pip install gdown\nthen re-run this script."
        )


def download_folder(url: str, dest: Path) -> None:
    """Download an entire Drive folder tree into `dest` using gdown."""
    dest.mkdir(parents=True, exist_ok=True)
    cmd = ["gdown", "--folder", url, "--output", str(dest), "--remaining-ok"]
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        sys.exit(
            "gdown failed. Common causes:\n"
            "  - the folder is not shared as 'anyone with the link'\n"
            "  - more than 50 files in one folder (gdown folder limit) — split it\n"
            "  - rate limiting (retry later)"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only",
        choices=sorted(SUBPROJECTS),
        help="fetch a single subproject instead of everything",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=REPO_ROOT,
        help="where to place the downloaded tree (default: repository root)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="show what would be fetched and exit",
    )
    args = parser.parse_args()

    targets = [args.only] if args.only else list(SUBPROJECTS)

    if args.list:
        print("Would fetch from the shared Drive folder:")
        for name in targets:
            print(f"  - {name}/")
        return

    check_gdown()
    print(f"Downloading {', '.join(targets)} into {args.dest} ...")
    # gdown fetches the whole folder; selective fetch requires the Drive folder
    # to expose per-subproject share links. With a single folder link we pull
    # the tree once and let the caller ignore what they do not need.
    download_folder(DRIVE_FOLDER_URL, args.dest)
    print("Done. Verify that files landed in the expected step directories.")


if __name__ == "__main__":
    main()
