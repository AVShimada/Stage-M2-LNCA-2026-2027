"""Comparer le MFM retenu et le contrôle sur cinq nouvelles graines du même scan."""
import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import commun as io
from audits.comparer_jumeau import render_comparison
from dynamique import TARGETS, fc, feature_set
from jumeau import setup, simulate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resultat", type=Path, help="Chemin du resultat.json du balayage retenu")
    args = parser.parse_args()
    result_file = args.resultat.resolve()
    sweep = result_file.parent.parent
    manifest = sweep / "manifest.json"
    state = io.read_json(manifest)
    if state["status"] != "complete":
        parser.error("Le balayage doit être terminé.")
    design = state["configuration"]
    cfg, model = design["pilote"], design["modele"]
    result = io.read_json(result_file)
    params = {label: {k: result[label][k] for k in ("G","sigma")}
              for label in ("controle","meilleur_sur_grille")}
    seeds = [201, 202, 203, 204, 205]
    if set(seeds).intersection([*model["fit_seeds"], *model["validation_seeds"]]):
        parser.error("Graines de confirmation déjà utilisées dans le balayage.")
    x, grid, n, sc_obs, names, files = io.load_inputs(cfg)
    sc, sc_files = setup(cfg, model, names)
    files.update(sc_files)
    for k, p in files.items():
        if io.digest(p) != state["inputs"][k]["sha256"]:
            raise ValueError(f"Entrée modifiée depuis le balayage : {k}")
    if cfg.get("fc_reference_z") and cfg["censor_mode"] == "none":
        files["fc_reference_z"] = io.source(cfg, "fc_reference_z")
        control_file = Path(design["controle_source"])
        if io.digest(control_file) != state["inputs"]["controle"]["sha256"]:
            raise ValueError("Manifeste du contrôle modifié depuis le balayage.")
        control_fc = io.read_json(control_file)["inputs"]["fc_reference_z"]
        if io.digest(files["fc_reference_z"]) != control_fc["sha256"]:
            raise ValueError("FC fournie différente de la référence du contrôle.")
        files["controle"] = control_file
    targets_file = result_file.parent / "cibles_empiriques.json"
    targets = io.read_json(targets_file)
    empirical = np.array([targets[k] for k in TARGETS])
    if not np.allclose(feature_set(x, grid, cfg["tr_seconds"])[0], empirical):
        raise ValueError("Les cibles empiriques recalculées diffèrent du balayage.")
    obs_edges = fc(x)[io.EDGE]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = io.stage_dir("12_evaluer_jumeaux") / io.run_dir(cfg, "", stamp).name
    files.update(balayage=manifest, resultat=result_file, cibles=targets_file, confirmation=__file__,
                 rendu=Path(__file__).with_name("comparer_jumeau.py"))
    settings = dict(pilote=cfg, modele_du_controle=model, parametres=params,
                    graines=seeds, graine_figures=seeds[0])
    print(out, flush=True)
    with io.record(out, "Confirmation du MFM recalibré sur le pilote", settings, files) as data:
        predictions = {}
        with (data/"confirmation.csv").open("x", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["modele","G","sigma","graine","erreur_10_cibles","r_entre_liens_FC",*TARGETS])
            for label, par in params.items():
                errors = []
                for seed in seeds:
                    y = simulate(sc, cfg, model, par, n, seed)
                    if label == "meilleur_sur_grille" and seed == seeds[0]:
                        np.save(data/"bold_simule.npy", y, allow_pickle=False)
                        _, q = render_comparison(out, cfg, dict(model, control_parameters=par),
                                                 x, y[grid], grid, sc_obs, sc)
                    else:
                        q = feature_set(y[grid], grid, cfg["tr_seconds"])[0]
                    edges = fc(y[grid])[io.EDGE]
                    r = np.corrcoef(edges, obs_edges)[0, 1] if min(np.ptp(edges), np.ptp(obs_edges)) > 0 else ""
                    loss = float(np.linalg.norm(q-empirical))
                    errors.append(loss)
                    writer.writerow([label, par["G"], par["sigma"], seed, loss, r, *q])
                    f.flush()
                    print(f"{label} · graine {seed} · erreur {loss:.4f}", file=sys.__stdout__, flush=True)
                predictions[label] = errors
        means = {label: float(np.mean(errors)) for label, errors in predictions.items()}
        wins = sum(b < a for a,b in zip(predictions["controle"], predictions["meilleur_sur_grille"]))
        summary = dict(erreur_moyenne=means, graines=seeds, n_graines_ameliorees=wins,
                       parametres=params, graine_figures=seeds[0])
        io.write_json(data/"confirmation.json", summary)
    print(summary)
    print(out/"figures")


if __name__ == "__main__":
    main()
