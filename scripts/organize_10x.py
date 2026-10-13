"""Group GSE121638 files by sample ID; preview by default and never overwrite."""

import argparse
import os
from pathlib import Path
import re


def organize(directory, apply=False):
    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    moves = []
    for source in sorted(root.iterdir()):
        if not source.is_file() or source.is_symlink():
            continue
        match = re.search(r"(?<![A-Za-z0-9])GU[0-9]+_[PT](?![A-Za-z0-9])", source.name)
        if match:
            target = root / match.group() / source.name
            if target.parent.is_symlink() or (target.parent.exists() and not target.parent.is_dir()):
                raise ValueError(f"Invalid sample directory: {target.parent}")
            if target.exists() or target.is_symlink():
                raise FileExistsError(f"Refusing to overwrite: {target}")
            moves.append((source, target))
    for source, target in moves:
        print(f"{source.name} -> {target.relative_to(root)}")
    if apply:
        for source, target in moves:
            target.parent.mkdir(exist_ok=True)
            # Link creation refuses an existing target, including races after preview.
            os.link(source, target)
            source.unlink()
    return moves


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--apply", action="store_true", help="Apply the previewed moves")
    args = parser.parse_args()
    organize(args.directory, args.apply)
