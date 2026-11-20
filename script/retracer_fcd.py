"""Tracer les deux échelles de référence à partir de la FCD à 60 s déjà calculée."""
import argparse
from pathlib import Path

import numpy as np

from commun import heatmap, read_json


def render(folder):
    state = read_json(folder / "manifest.json")
    if state["status"] != "complete":
        raise ValueError("Le lancement doit être terminé pour tracer sa FCD.")
    config = state["configuration"]
    fcd = np.load(folder / "data/fcd.npy", allow_pickle=False)
    if fcd.ndim != 2 or fcd.shape[0] != fcd.shape[1] or not np.isfinite(fcd).all():
        raise ValueError("Une matrice de FCD carrée et finie est attendue.")
    title = f"FCD · {config['subject']} / {config['session']} / {config['version']} · 60 s"
    for name, label in (("fcd_m1", "échelle M1"), ("fcd_article", "échelle article")):
        heatmap(folder, name, fcd, title+" · "+label, "Fenêtre", style=name)
        print(folder / "figures" / (name+".png"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lancement", type=Path, help="Dossier contenant manifest.json et data/fcd.npy")
    render(parser.parse_args().lancement)
