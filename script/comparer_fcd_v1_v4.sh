#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1

"${PYTHON_BIN:-$HOME/.venvs/lnca-1000brains/bin/python}" - "${1:-config/pilote.json}" <<'PY'
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

# Donner la priorité aux modules installés dans script/ pour ce projet.
sys.path.insert(0, str((Path.cwd() / "script").resolve()))
from commun import array_file, heatmap, read_json, record, run_dir, source, write_json
from dynamique import unit, windows

cfg = read_json(sys.argv[1])
bold_v4 = source(cfg, "bold")
suffix = "_V4.txt"
if cfg["version"] != "V4" or not bold_v4.name.endswith(suffix):
    raise ValueError("Utiliser la configuration du pilote et sa série BOLD V4.")
if cfg.get("censor_mode") != "none":
    raise ValueError("Ce lanceur utilise des volumes consécutifs. Fournir une grille commune si des volumes ont été retirés.")

bold_v1 = bold_v4.with_name(bold_v4.name.removesuffix(suffix) + "_V1.txt")
inputs = {"V1": bold_v1, "V4": bold_v4}
signals = {version: array_file(path) for version, path in inputs.items()}
if cfg["bold_orientation"] == "regions_time":
    signals = {version: signal.T for version, signal in signals.items()}
elif cfg["bold_orientation"] != "time_regions":
    raise ValueError("Renseigner l’orientation des séries BOLD.")
shapes = {version: signal.shape for version, signal in signals.items()}
if shapes["V1"] != shapes["V4"]:
    raise ValueError(f"Dimensions V1 et V4 différentes : {shapes}. Vérifier les indices des volumes d’origine.")

width, step = 10, 1
signal = signals["V4"]
if signal.ndim != 2 or signal.shape[1] != 100 or len(signal) < 2 * width:
    raise ValueError(f"Séries attendues : au moins {2 * width} volumes × 100 régions, reçu {signal.shape}.")
if not all(np.isfinite(signal).all() for signal in signals.values()):
    raise ValueError("Valeurs non finies dans les séries BOLD.")
grid = np.arange(len(signal))
limits = (0, 0.5)
run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
out = run_dir(cfg, "02_comparaison_fcd_v1_v4", run)
settings = {
    "width_TR": width,
    "step_TR": step,
    "n_volumes_per_version": len(grid),
    "display_limits": limits,
    "hypothese": "Les lignes V1 et V4 correspondent aux mêmes volumes consécutifs du scan.",
}
with record(out, "Comparaison des FCD V1 et V4", cfg, inputs, extra=settings) as data:
    summary = {}
    for version, signal in signals.items():
        stream, starts = windows(signal, grid, width, step)
        normalized = unit(stream, axis=1)
        fcd = np.clip(normalized @ normalized.T, -1, 1)
        first, second = np.triu_indices(len(starts), k=1)
        nonoverlapping = starts[second] - starts[first] >= width
        chi = fcd[first[nonoverlapping], second[nonoverlapping]]
        summary[version] = {
            "n_fenetres": len(starts),
            "n_paires_sans_chevauchement": len(chi),
            "chi_mediane_sans_chevauchement": float(np.median(chi)),
        }
        np.save(data / f"fcd_{version.lower()}.npy", fcd)
        title = f"FCD {version} · {cfg['subject']} / {cfg['session']} · 10 volumes"
        heatmap(out, f"fcd_{version.lower()}", fcd, title, "Fenêtre",
                style="fcd_article", limits=limits)
    np.save(data / "starts_TR.npy", starts)
    write_json(data / "resume.json", summary)
print(f"V1 et V4 : {len(grid)} volumes, {len(starts)} fenêtres par version.")
for version, stats in summary.items():
    print(f"{version} : χ médian sans chevauchement = {stats['chi_mediane_sans_chevauchement']:.4f}")
print(out)
PY
