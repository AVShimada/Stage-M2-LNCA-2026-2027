"""Retracer la MC du pilote avec les liens et les bornes de la carte non triée."""
import argparse
import csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import commun as io


def render(out, mc, labels):
    labels = np.asarray(labels)
    if mc.shape != (len(labels), len(labels)):
        raise ValueError("Les dimensions de la MC et le nombre de liens classés diffèrent.")
    view = mc[::5, ::5]
    if not np.isfinite(view).all():
        raise ValueError("La MC affichée contient des valeurs non finies.")
    vmax = max(np.finfo(float).eps, float(np.percentile(view, 99)))
    groups = labels[::5]
    order = np.argsort(groups, kind="stable")
    sorted_view = view[np.ix_(order, order)]
    boundaries = np.flatnonzero(np.diff(groups[order])) + 0.5

    def draw(fig, ax):
        im = ax.imshow(sorted_view, origin="lower", aspect="auto", interpolation="nearest",
                       cmap="turbo", vmin=0, vmax=vmax)
        for edge in boundaries:
            ax.axhline(edge, linestyle=":", linewidth=0.7, color="white")
            ax.axvline(edge, linestyle=":", linewidth=0.7, color="white")
        ax.set(title="Méta-connectivité triée par module", xlabel="Rang du lien affiché",
               ylabel="Rang du lien affiché")
        tails = (view.min() < 0, view.max() > vmax)
        extend = {(False, False): "neither", (True, False): "min",
                  (False, True): "max", (True, True): "both"}[tails]
        fig.colorbar(im, ax=ax, extend=extend, label="Corrélation de méta-connectivité")

    # Finir les deux rendus avant de remplacer les figures existantes.
    dest = out / "figures"
    with TemporaryDirectory(dir=dest) as tmp:
        staging = Path(tmp)
        (staging / "figures").mkdir()
        io.figure(staging, "modules_mc", draw, dark=True)
        for ext in ("png", "pdf"):
            file = f"modules_mc.{ext}"
            (staging / "figures" / file).replace(dest / file)
    return len(view), vmax


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, help="Identifiant UTC du pilote sub-0001 / ses-1 / V4")
    args = parser.parse_args()
    tag = f"{io.token(args.run)}_sub-0001_ses-1_V4"
    source = io.stage_dir("02_dynamique_m1") / tag
    out = io.stage_dir("07_modules_mc") / tag
    for folder in (source, out):
        if io.read_json(folder / "manifest.json")["status"] != "complete":
            parser.error(f"Étape encore incomplète : {folder}")
    mc = np.load(source / "data/mc_20s.npy", mmap_mode="r", allow_pickle=False)
    with (out / "data/modules.csv").open(encoding="utf-8-sig", newline="") as f:
        pairs = sorted((int(r["edge_0"]), int(r["module"])) for r in csv.DictReader(f))
    if [i for i, _ in pairs] != list(range(len(mc))):
        parser.error("modules.csv doit contenir chaque indice de lien exactement une fois.")
    n, high = render(out, mc, [label for _, label in pairs])
    print(f"{n} liens affichés, échelle commune de 0 à {high:.6g}.")
    print(out / "figures/modules_mc.png")
    print(out / "figures/modules_mc.pdf")


if __name__ == "__main__":
    main()
