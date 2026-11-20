#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=script OPENBLAS_NUM_THREADS=1

"${PYTHON_BIN:-$HOME/.venvs/lnca-1000brains/bin/python}" - "${1:-config/pilote.json}" <<'PY'
from datetime import datetime, timezone
import sys

import numpy as np

from commun import EDGE, array_file, figure, read_json, record, run_dir, source

# Minima et maxima communs aux deux groupes de la figure M1 non triée.
FC_LIMITS_M1 = (-0.40065664835678494, 1.2926806413667378)
SC_LIMITS_M1 = (-0.08296555948553058, 4.607515607717042)


def read_matrix(path):
    matrix = array_file(path)
    if (matrix.shape != (100, 100) or not np.isfinite(matrix[EDGE]).all()
            or not np.allclose(matrix[EDGE], matrix.T[EDGE], atol=1e-8)):
        raise ValueError(f"Matrice 100 × 100 symétrique et finie hors diagonale attendue : {path}")
    return matrix


def plot_matrix(out, name, matrix, title, limits, show_tails):
    display = matrix.copy()
    np.fill_diagonal(display, 0)

    def draw(fig, ax):
        image = ax.imshow(display, cmap="turbo", origin="upper", aspect="equal",
                          interpolation="nearest", vmin=limits[0], vmax=limits[1])
        ax.set(title=title, xlabel="Nœuds", ylabel="Nœuds")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, extend=show_tails)

    figure(out, name, draw)


config = read_json(sys.argv[1])
version4 = source(config, "fc_reference_z")
suffix = "_V4.txt"
if config["version"] != "V4" or not version4.name.endswith(suffix):
    raise ValueError("Ce contrôle attend la configuration V4 et son fichier FC Fisher z.")
version1 = version4.with_name(version4.name.removesuffix(suffix) + "_V1.txt")
inputs = {"sc": source(config, "sc"), "fc_v1": version1, "fc_v4": version4}
matrices = {name: read_matrix(path) for name, path in inputs.items()}
p2, p98 = map(float, np.percentile(matrices["fc_v4"][EDGE], [2, 98]))
v4_limits = (p2, p98) if p2 < p98 else FC_LIMITS_M1
plots = {
    "sc": ("sc", f"SC · {config['sc_weights']}", SC_LIMITS_M1),
    "fc_v1": ("fc_v1", "FC V1 · z de Fisher", FC_LIMITS_M1),
    "fc_v4": ("fc_v4", "FC V4 · z de Fisher", FC_LIMITS_M1),
    "fc_v4_adaptee": ("fc_v4", "FC V4 · z de Fisher · échelle adaptée", v4_limits),
}
tails = {}
for name, (matrix_name, _, (low, high)) in plots.items():
    links = matrices[matrix_name][EDGE]
    tails[name] = {
        "n_links": len(links),
        "n_below": int(np.count_nonzero(links < low)),
        "n_above": int(np.count_nonzero(links > high)),
    }

run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
out = run_dir(config, "01_apercu_sc_fc", run)
scan = f"{config['subject']} / {config['session']}"
reference = {"reference": "group_average_matrices_FC_SC.png (M1)",
             "fc_limits": FC_LIMITS_M1, "sc_limits": SC_LIMITS_M1,
             "fc_v4_adaptee_limits": v4_limits, "fc_v4_adaptee_percentiles": [2, 98],
             "links_outside_display_limits": tails}
with record(out, "Aperçu SC et FC", config, inputs, extra=reference) as data:
    for name, matrix in matrices.items():
        np.save(data / f"{name}.npy", matrix)
    extend = {(False, False): "neither", (True, False): "min",
              (False, True): "max", (True, True): "both"}
    for name, (matrix_name, label, limits) in plots.items():
        counts = tails[name]
        show_tails = extend[bool(counts["n_below"]), bool(counts["n_above"])]
        plot_matrix(out, name, matrices[matrix_name], f"{label} · {scan}", limits, show_tails)
print(out / "figures")
PY
