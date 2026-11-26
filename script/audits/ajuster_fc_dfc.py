"""Ajuster la FC du pilote sous contrainte sur les dix cibles de dFC.

La référence vient du manifeste de confirmation. La grille varie G et w à
sigma fixé. La RMSE compare les 4 950 liens de Pearson calculés sur les BOLD.
Les dix quantiles gardent la définition de feature_set. La marge relative
s'applique à l'erreur moyenne de la référence, séparément en ajustement et
en validation. Les figures illustrent le réglage retenu pour la graine 301.
"""
import argparse
import csv
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import commun as io
from audits.comparer_jumeau import PILOTE, render_comparison
from dynamique import TARGETS, fc, feature_set
from jumeau import setup, simulate


METRICS = ("rmse_fc", "erreur_10_cibles", "mediane_fc_r", "r_entre_liens_fc")
VALIDATION = (301, 302, 303, 304, 305)


def compare_series(series, grid, tr, empirical_edges, empirical_targets):
    edges = fc(series)[io.EDGE]
    targets = feature_set(series, grid, tr)[0]
    rmse = np.sqrt(np.mean((edges - empirical_edges) ** 2))
    distance = np.linalg.norm(targets - empirical_targets)
    correlation = None
    if np.ptp(edges) > 0 and np.ptp(empirical_edges) > 0:
        correlation = float(np.corrcoef(edges, empirical_edges)[0, 1])
    return dict(rmse_fc=float(rmse), erreur_10_cibles=float(distance),
                mediane_fc_r=float(np.median(edges)), r_entre_liens_fc=correlation,
                **dict(zip(TARGETS, targets.tolist())))


def select_candidate(scores, reference, margin):
    limit = scores[reference]["erreur_10_cibles"] * (1 + margin)
    eligible = [key for key, v in scores.items() if v and v["erreur_10_cibles"] <= limit]
    # À RMSE égale à celle de la référence, conserver la référence.
    return min(eligible, key=lambda k: (scores[k]["rmse_fc"], k != reference,
                                        scores[k]["erreur_10_cibles"], k))


def plot_grid(out, training, reference, best, margin):
    def draw(fig, axis):
        axis.remove()
        axes = fig.subplots(1, 2)
        for ax, metric, label in zip(axes, METRICS[:2], ("RMSE des liens FC (r)", "Erreur sur les dix cibles dFC")):
            for w in sorted({pair[1] for pair in training}):
                keys = sorted(pair for pair in training if pair[1] == w)
                ax.plot([p[0] for p in keys], [training[p].get(metric, np.nan) for p in keys],
                        "o-", label=f"w = {w:g}")
            ax.scatter(best[0], training[best][metric], s=90, marker="s", facecolors="none",
                       edgecolors="black", zorder=4, label="Candidat retenu en ajustement")
            ax.set(xlabel="Couplage global G", ylabel=label)
        axes[1].axhline(training[reference]["erreur_10_cibles"] * (1 + margin),
                        color="black", linestyle=":", label="Limite dFC")
        axes[1].legend(fontsize=8)
        fig.suptitle("FC et dFC du pilote · moyennes sur les graines d’ajustement")
    io.figure(out, "06_balayage_fc_dfc", draw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("confirmation", type=Path, help="Dossier du run de confirmation de la dFC")
    parser.add_argument("--g", nargs="+", type=float, default=[0.60, 0.75, 0.90])
    parser.add_argument("--w", nargs="+", type=float, default=[0.90, 1.00, 1.10])
    parser.add_argument("--marge-dfc", type=float, default=0.10,
                        help="Hausse relative maximale de l'erreur dynamique (défaut 0.10)")
    args = parser.parse_args()
    manifest = args.confirmation.resolve() / "manifest.json"
    state = io.read_json(manifest)
    design = state["configuration"]
    cfg, model = design["pilote"], design["modele_du_controle"]
    if state["status"] != "complete" or tuple(cfg[k] for k in ("subject", "session", "version")) != PILOTE:
        parser.error("Fournir le dossier de confirmation terminé du pilote V4.")
    confirmed = design["parametres"]["meilleur_sur_grille"]
    sigma = confirmed["sigma"]
    reference = (confirmed["G"], model["fixed"]["w"])
    candidates = [reference, *sorted(set(product(args.g, args.w)) - {reference})]
    if not all(np.isfinite(v) and v >= 0 for v in [*args.g, *args.w, args.marge_dfc]):
        parser.error("G, w et la marge doivent être finis et positifs ou nuls.")
    fit = model["fit_seeds"]
    used_seeds = set(fit) | set(model["validation_seeds"]) | set(design["graines"])
    if not fit or len(set(fit)) != len(fit) or used_seeds.intersection(VALIDATION):
        parser.error("Vérifier les graines d'ajustement et l'indépendance de la nouvelle validation.")

    x, grid, n, sc_obs, names, files = io.load_inputs(cfg)
    sc, sc_files = setup(cfg, model, names)
    files.update(sc_files)
    if cfg.get("fc_reference_z") and cfg["censor_mode"] == "none":
        files["fc_reference_z"] = io.source(cfg, "fc_reference_z")
    files["cibles"] = Path(state["inputs"]["cibles"]["path"])
    for key, path in files.items():
        if io.digest(path) != state["inputs"][key]["sha256"]:
            raise ValueError(f"Entrée différente de la confirmation : {key}")
    target = feature_set(x, grid, cfg["tr_seconds"])[0]
    saved = io.read_json(files["cibles"])
    if not np.allclose(target, [saved[k] for k in TARGETS]):
        raise ValueError("Les cibles dynamiques du pilote ont changé.")
    empirical_edges = fc(x)[io.EDGE]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = io.stage_dir("10_ajuster_mfm") / "fc_dfc" / io.run_dir(cfg, "", stamp).name
    settings = dict(pilote=cfg, modele=model, reference=dict(G=reference[0], w=reference[1], sigma=sigma),
                    grille_G_w=candidates, graines_ajustement=fit, graines_validation=VALIDATION,
                    marge_relative_dfc=args.marge_dfc, graine_figures=VALIDATION[0])
    files.update(manifeste_confirmation=manifest, script_ajustement=__file__,
                 rendu=Path(__file__).with_name("comparer_jumeau.py"))
    print(out, flush=True)
    with io.record(out, "MFM du pilote · ajustement conjoint FC et dFC", settings, files) as data:
        examples = {}
        averages = []
        with (data / "simulations.csv").open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, ["phase", "G", "w", "sigma", "graine", "statut", *METRICS, *TARGETS])
            writer.writeheader()

            def evaluate(pair, seeds, phase):
                params = dict(G=pair[0], w=pair[1], sigma=sigma)
                rows = []
                for seed in seeds:
                    info = dict(phase=phase, **params, graine=seed)
                    try:
                        y = simulate(sc, cfg, model, params, n, seed)[grid]
                        row = compare_series(y, grid, cfg["tr_seconds"], empirical_edges, target)
                    except ValueError as exc:
                        if not str(exc).startswith(("Intégration instable", "Corrélation indéfinie")):
                            raise
                        writer.writerow(dict(**info, statut=str(exc)))
                    else:
                        rows.append(row)
                        writer.writerow(dict(**info, statut="ok", **row))
                        if phase == "validation" and seed == VALIDATION[0]:
                            examples[pair] = y
                    stream.flush()
                mean = {}
                if len(rows) == len(seeds):
                    for metric in METRICS:
                        values = [r[metric] for r in rows]
                        mean[metric] = float(np.mean(values)) if None not in values else None
                averages.append([phase, *pair, sigma, len(rows), len(seeds), *[mean.get(k) for k in METRICS]])
                message = (f"{phase} G={pair[0]:g}, w={pair[1]:g} · " +
                           (f"RMSE FC {mean['rmse_fc']:.4f} · erreur dFC {mean['erreur_10_cibles']:.4f}"
                            if mean else f"{len(rows)}/{len(seeds)} simulations exploitables"))
                print(message, flush=True)
                print(message, file=sys.__stdout__, flush=True)
                return mean

            training = {pair: evaluate(pair, fit, "ajustement") for pair in candidates}
            if not training[reference]:
                raise ValueError("Référence non estimable sur toutes les graines d'ajustement.")
            best = select_candidate(training, reference, args.marge_dfc)
            validation = {pair: evaluate(pair, VALIDATION, "validation")
                          for pair in dict.fromkeys((reference, best))}
        io.write_table(data / "grille.csv", ["phase", "G", "w", "sigma", "n_ok", "n_prevu", *METRICS], averages)
        if not validation[reference]:
            raise ValueError("Référence non estimable sur les cinq graines de validation.")
        selected = select_candidate(validation, reference, args.marge_dfc)

        def report(pair):
            return dict(G=pair[0], w=pair[1], sigma=sigma, ajustement=training[pair],
                        validation=validation[pair] or None)

        result = dict(reference=report(reference), candidat_ajustement=report(best), retenu=report(selected),
                      changement_confirme=selected != reference, marge_relative_dfc=args.marge_dfc,
                      fc_empirique_mediane_r=float(np.median(empirical_edges)), graine_figures=VALIDATION[0])
        np.save(data / "bold_simule_volumes_retenus.npy", examples[selected], allow_pickle=False)
        plot_grid(out, training, reference, best, args.marge_dfc)
        render_comparison(out, cfg, model, x, examples[selected], grid, sc_obs, sc)
        io.write_json(data / "resultat.json", result)
    print("Réglage retenu", result["retenu"], flush=True)
    print(out / "data" / "resultat.json")


if __name__ == "__main__":
    main()
