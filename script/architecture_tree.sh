#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
"${PYTHON_BIN:-python3}" - <<'PY'
from pathlib import Path

root = Path.cwd()
report = root / "architecture_projet.txt"
omit = {"1kBrains", ".git", "__pycache__", ".venv", "tmp"}
lines = [str(root)]

def visit(folder, prefix="", depth=1):
    entries = []
    for p in folder.iterdir():
        if p.name in omit or p.name.startswith("._") or p.name == ".DS_Store" or p.suffix in {".gz", ".fna", ".faa"}:
            continue
        if folder == root and p.is_dir() and p.name != "output":
            if p.name.startswith(("output_", "output-", "output.")) or (p / "00_controle").is_dir() or (p / "02_pilote").is_dir():
                continue
        entries.append(p)
    for index, p in enumerate(sorted(entries, key=lambda item: item.name.casefold())):
        last = index == len(entries) - 1
        lines.append(prefix + ("└── " if last else "├── ") + p.name)
        if p.is_dir() and not p.is_symlink() and p.name != "archives" and depth < 5:
            visit(p, prefix + ("    " if last else "│   "), depth + 1)

visit(root)
report.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(report)
PY
