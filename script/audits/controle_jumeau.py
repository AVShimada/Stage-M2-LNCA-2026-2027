"""Simuler une fois le MFM du pilote V4 et comparer ses dix cibles empiriques."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
from shutil import copyfile
from time import perf_counter

import numpy as np

import commun as io
from audits.comparer_jumeau import render_comparison
from dynamique import TARGETS, feature_set, volumes
from jumeau import setup, simulate


def empirical_run(requested):
    base = io.stage_dir("03_cibles_2026")
    suffix = "_sub-0001_ses-1_V4"
    paths = [base / (io.token(requested) + suffix) / "manifest.json"] if requested else base.glob(f"*{suffix}/manifest.json")
    completed = []
    for path in paths:
        if path.is_file():
            state = io.read_json(path)
            cfg = state["configuration"]
            if state["status"] == "complete" and (cfg["subject"], cfg["session"], cfg["version"]) == ("sub-0001", "ses-1", "V4"):
                completed.append((state["started"], path, state))
    if not completed:
        raise ValueError(f"Aucun calcul de cibles terminé pour sub-0001 / ses-1 / V4 dans {base}.")
    _, path, state = max(completed)
    return path, state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", help="Identifiant du calcul des cibles, sans le suffixe sujet/session/version")
    parser.add_argument("--run", help="Identifiant du lancement de sortie")
    parser.add_argument("--simulation-seule", action="store_true",
                        help="Réserver le rendu des figures à l’étape de comparaison")
    args = parser.parse_args()

    manifest_path, state = empirical_run(args.source_run)
    cfg = state["configuration"]
    x, grid, n, sc_observed, names, inputs = io.load_inputs(cfg)
    for key, path in inputs.items():
        if state["inputs"].get(key, {}).get("sha256") != io.digest(path):
            raise ValueError(f"Entrée différente du calcul des cibles : {key}")
    targets_path = manifest_path.parent / "data/targets.json"
    targets = io.read_json(targets_path)
    observed = np.array([targets[key] for key in TARGETS])
    check, *_ = feature_set(x, grid, cfg["tr_seconds"])
    if not np.allclose(observed, check, rtol=1e-8, atol=1e-10):
        raise ValueError("Les cibles enregistrées diffèrent du calcul sur le BOLD du pilote.")

    template = io.ROOT / "config/jumeau.json"
    model = io.read_json(template)
    defaults = {"sc_model_labels": cfg["labels_sc"], "sc_normalization": "column_max",
                "neural_noise_convention": "diffusion_sqrt_dt", "global_signal_regression": False,
                "observation_band_hz": []}
    applied = {key: value for key, value in defaults.items() if model.get(key) in (None, "")}
    model.update(applied)
    parameters = {key: model["control_parameters"].get(key) for key in ("G", "sigma")}
    for key, value in (("G", 0.05), ("sigma", 0.03)):
        if parameters[key] is None:
            parameters[key] = value
            applied[key] = value
    model["control_parameters"] = parameters
    sc, model_inputs = setup(cfg, model, names)
    seed = model["fit_seeds"][0]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = io.run_dir(cfg, "09_jumeau_controle_technique", args.run or stamp)
    settings = {"pilote": cfg, "modele": model, "valeurs_de_controle": applied, "graine": seed}
    files = {**inputs, **model_inputs, "modele_source": template, "cibles": targets_path,
             "manifest_cibles": manifest_path, "script_controle": __file__,
             "script_comparaison": Path(__file__).with_name("comparer_jumeau.py")}
    if cfg.get("fc_reference_z") and cfg["censor_mode"] == "none":
        files["fc_reference_z"] = io.source(cfg, "fc_reference_z")
        recorded_fc = state["inputs"].get("fc_reference_z")
        if recorded_fc and io.digest(files["fc_reference_z"]) != recorded_fc["sha256"]:
            raise ValueError("La FC fournie a changé depuis le calcul des cibles.")

    print(f"Cibles du run {manifest_path.parent.name}. Simulation MFM avec {parameters}.", flush=True)
    with io.record(out, "Contrôle technique du jumeau V4", settings, files) as data:
        start = perf_counter()
        y = simulate(sc, cfg, model, parameters, n, seed)
        if y.shape != (n, len(names)) or not np.isfinite(y).all():
            raise ValueError("Le BOLD simulé comporte une dimension ou une valeur inattendue.")
        elapsed = perf_counter() - start
        np.save(data / "bold_simule.npy", y, allow_pickle=False)
        reference_fcd = "fcd_simulee_cibles_60s.npy"
        if args.simulation_seule:
            simulated, _, _, fcd, _, speeds, _ = feature_set(y[grid], grid, cfg["tr_seconds"])
            np.save(data / reference_fcd, fcd, allow_pickle=False)
            comparison = {"fcd_figure_shape": None, "fenetre_fcd_figure_volumes": None, "n_vitesses_mfm": len(speeds)}
        else:
            _, simulated = render_comparison(out, cfg, model, x, y[grid], grid, sc_observed, sc)
            comparison = io.read_json(data / "comparaison.json")
        copyfile(data / reference_fcd, data / "fcd_simulee.npy")
        target_shape = list(np.load(data / reference_fcd, mmap_mode="r", allow_pickle=False).shape)
        figure_shape = comparison["fcd_figure_shape"]
        width = volumes(60, cfg["tr_seconds"])
        io.write_json(data / "modele_controle.json", model)
        io.write_json(data / "controle.json", {
            "bold_shape": list(y.shape),
            "fcd_shape": target_shape,
            "fcd_cibles_shape": target_shape,
            "fenetre_cibles_volumes": width,
            "fichier_fcd_cibles_simulee": reference_fcd,
            "alias_fichiers": {"fcd_simulee.npy": reference_fcd},
            "fcd_figure_shape": figure_shape,
            "fenetre_fcd_figure_volumes": comparison["fenetre_fcd_figure_volumes"],
            "fichiers_fcd_figures": ["fcd_empirique_10vol.npy", "fcd_simulee_10vol.npy"] if figure_shape else [],
            "n_vitesses": comparison["n_vitesses_mfm"], "secondes_calcul": elapsed,
            "distance_euclidienne_cibles": float(np.linalg.norm(simulated - observed)),
        })
        # Conserver les en-têtes des précédents contrôles techniques.
        io.write_table(data / "comparaison_cibles.csv", ["cible", "empirique", "simulee", "ecart"],
                       zip(TARGETS, observed, simulated, simulated - observed))
    print(f"BOLD simulé {y.shape} | {elapsed:.1f} s de simulation")
    if figure_shape:
        print(f"FCD affichée · 10 volumes · {tuple(figure_shape)}")
    print(f"FCD des cibles · fenêtre de référence 60 s ({width} volumes) · {tuple(target_shape)}")
    print(out)


if __name__ == "__main__":
    main()
