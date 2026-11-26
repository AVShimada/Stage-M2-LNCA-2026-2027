"""Balayage de G et sigma sur le MFM d'un pilote déjà simulé."""
import argparse
import csv
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import commun as io
from audits.comparer_jumeau import control_manifest
from dynamique import TARGETS, fc, feature_set
from jumeau import setup, simulate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", required=True, help="Identifiant du contrôle technique")
    parser.add_argument("--g", nargs="+", type=float, default=[0.025, 0.05, 0.1, 0.2, 0.4])
    parser.add_argument("--sigma", nargs="+", type=float, default=[0.015, 0.03, 0.06])
    args = parser.parse_args()
    path, state, cfg, model, _, _ = control_manifest(args.source_run)
    x, grid, n, _, names, files = io.load_inputs(cfg)
    sc, model_files = setup(cfg, model, names)
    files.update(model_files)
    for key, p in files.items():
        if io.digest(p) != state["inputs"][key]["sha256"]:
            raise ValueError(f"Entrée différente du contrôle : {key}")
    reference = tuple(float(model["control_parameters"][k]) for k in ("G", "sigma"))
    gs = sorted(set([reference[0], *args.g]))
    sigmas = sorted(set([reference[1], *args.sigma]))
    if not all(np.isfinite(v) and v >= 0 for v in gs + sigmas):
        parser.error("Les paramètres doivent être finis et positifs ou nuls.")
    fit, valid = model["fit_seeds"], model["validation_seeds"]
    if not fit or not valid or set(fit) & set(valid) or len(set(fit)) != len(fit) or len(set(valid)) != len(valid):
        parser.error("Les graines doivent être distinctes dans et entre les deux séries.")
    target = feature_set(x, grid, cfg["tr_seconds"])[0]
    emp_fc = fc(x)[io.EDGE]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = io.stage_dir("10_ajuster_mfm") / "balayages" / io.run_dir(cfg, "", stamp).name
    design = dict(pilote=cfg, modele=model, G=gs, sigma=sigmas, controle_source=str(path))
    files.update(controle=path, script_balayage=__file__)
    print(out, flush=True)
    with io.record(out, "Balayage du MFM et validation sur graines indépendantes", design, files) as data:
        io.write_json(data/"cibles_empiriques.json", dict(zip(TARGETS, target.tolist())))
        with (data/"simulations.csv").open("x", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["phase","G","sigma","graine","statut","erreur_10_cibles",
                             "correlation_liens_FC_pearson","FC_pearson_mediane",*TARGETS])

            def score(pair, seeds, phase):
                g, sigma = pair
                errors = []
                for seed in seeds:
                    prefix = [phase, g, sigma, seed]
                    try:
                        y = simulate(sc, cfg, model, {"G":g, "sigma":sigma}, n, seed)[grid]
                        pred = feature_set(y, grid, cfg["tr_seconds"])[0]
                        edges = fc(y)[io.EDGE]
                    except ValueError as exc:
                        if not str(exc).startswith(("Intégration instable","Corrélation indéfinie")):
                            raise
                        writer.writerow([*prefix,str(exc),*([""]*(3+len(TARGETS)))])
                    else:
                        error = float(np.linalg.norm(pred-target))
                        correlation = np.corrcoef(edges, emp_fc)[0,1] if min(np.ptp(edges), np.ptp(emp_fc)) > 0 else ""
                        writer.writerow([*prefix,"ok",error,correlation,float(np.median(edges)),*pred])
                        errors.append(error)
                    f.flush()
                mean = float(np.mean(errors)) if len(errors)==len(seeds) else None
                msg = f"{phase} G={g:g}, sigma={sigma:g} : {mean if mean is not None else 'non estimable'}"
                print(msg, flush=True)
                print(msg, file=sys.__stdout__, flush=True)
                return mean

            scores = {pair: score(pair, fit, "ajustement") for pair in product(gs, sigmas)}
            finite = [pair for pair, error in scores.items() if error is not None and np.isfinite(error)]
            if not finite:
                raise ValueError("Aucun point de grille estimable sur toutes les graines.")
            best = min(finite, key=scores.get)
            validation = {p:score(p,valid,"validation") for p in dict.fromkeys([reference,best])}
        def report(p):
            return {"G":p[0],"sigma":p[1],"erreur_ajustement":scores[p],"erreur_validation":validation[p]}
        result = {"controle":report(reference),"meilleur_sur_grille":report(best)}
        io.write_json(data/"resultat.json", result)
        def draw(fig, ax):
            matrix = np.array([[scores[g,s] if scores[g,s] is not None else np.nan for g in gs] for s in sigmas])
            im = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis_r", interpolation="nearest")
            ax.set_xticks(range(len(gs)), [f"{v:g}" for v in gs])
            ax.set_yticks(range(len(sigmas)), [f"{v:g}" for v in sigmas])
            ax.set(xlabel="Couplage global G", ylabel="Amplitude du bruit σ", title="Erreur moyenne sur les dix cibles de dFC")
            fig.colorbar(im, ax=ax, label="Distance euclidienne moyenne")
        io.figure(out, "erreur_G_sigma", draw)
    print(result)
    print(out/"data/resultat.json")


if __name__ == "__main__":
    main()
