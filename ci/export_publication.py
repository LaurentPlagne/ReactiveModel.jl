"""Create a reviewable source snapshot, excluding Git history and private audits.

    python3 ci/export_publication.py /absolute/new-or-empty-directory

The output must be outside the checkout. This copies the current files and does
not establish their provenance; review inputs and measured source bundles first.
It neither initializes a repository nor changes the original history.
"""
from pathlib import Path
import argparse
import hashlib
import shutil


# Historical evidence stays outside the current publication source tree.
EXCLUDED_ARTIFACTS = (
    "compilation-source-2026-11-13.tar.gz",
    "independent-source-2026-11-13.tar.gz",
    "compilation-results-2026-11-13.toml",
    "independent-results-2026-11-13.toml",
)


def export(destination):
    root = Path(__file__).resolve().parents[1]
    destination = destination.resolve()
    if destination == root or root in destination.parents:
        raise ValueError("choose an output directory outside the checkout")
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise ValueError("output must be a new or empty directory")
    destination.mkdir(parents=True, exist_ok=True)
    ignore = shutil.ignore_patterns(".git", "build", "__pycache__", "*.pyc", "*.jl.cov", "*.jl.mem", *EXCLUDED_ARTIFACTS)
    for name in (".github", "ci", "src", "ext", "test", "examples", "docs"):
        shutil.copytree(root / name, destination / name, ignore=ignore)
    for name in (".gitignore", "Project.toml", "Manifest.toml", "README.md", "CONTRIBUTING.md", "LICENSE"):
        shutil.copy2(root / name, destination / name)
    benchmarks = destination / "benchmarks"
    benchmarks.mkdir()
    for file in (root / "benchmarks").iterdir():
        if file.is_file() and not file.name.startswith("render_") and file.suffix in (".jl", ".py", ".md", ".toml"):
            shutil.copy2(file, benchmarks / file.name)
    # The documentation link checker is part of CI; other audit files are not
    # part of the publication snapshot.
    (destination / "audit").mkdir()
    shutil.copy2(root / "audit/check-documentation-links.py", destination / "audit/check-documentation-links.py")
    files = sorted(p for p in destination.rglob("*") if p.is_file())
    checksums = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(destination).as_posix()}\n" for p in files]
    (destination / "SHA256SUMS").write_text("".join(checksums))
    print(f"Exported {len(files)} files and SHA256SUMS to {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    try:
        export(args.destination)
    except ValueError as error:
        parser.error(str(error))
