"""Comparer les figures empiriques V4 et celles du MFM déjà simulé."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy.stats import gaussian_kde

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import commun as io
from dynamique import PERCENTILES, TARGETS, fc, feature_set, volumes
from jumeau import setup
from rendu_pilote import fc_display_limits, fcd_ten_volumes, plot_static_matrix


PILOTE = ("sub-0001", "ses-1", "V4")
EMP_COLOR, MFM_COLOR = "#222222", io.MODEL_COLORS["mfm"]


def control_manifest(run):
    suffix = "_" + "_".join(PILOTE)
    stages = ("09_jumeau_controle_technique", "09_simulation")
    paths = (path for stage in stages for path in io.stage_dir(stage).glob(f"*{suffix}/manifest.json"))
    runs = []
    for path in paths:
        if run and (path.parent.name not in (run, run + suffix)
                    and Path(run).resolve() != path.parent.resolve()):
            continue
        state = io.read_json(path)
        if state["status"] != "complete":
            continue
        if "model_configuration" in state:
            cfg, model = state["configuration"], state["model_configuration"]
            seed = model["fit_seeds"][0]
            bold = path.parent / "data/simulated_bold.npy"
        else:
            settings = state["configuration"]
            cfg, model, seed = settings["pilote"], settings["modele"], settings["graine"]
            bold = path.parent / "data/bold_simule.npy"
        if (cfg["subject"], cfg["session"], cfg["version"]) == PILOTE:
            runs.append((state["started"], path, state, cfg, model, seed, bold))
    if not runs:
        raise ValueError("Aucune simulation MFM terminée pour ce pilote sous output/03_jumeau/.")
    if run and len(runs) > 1:
        raise ValueError("Deux simulations portent ce nom. Donner le chemin complet de la simulation avec --source-run.")
    return max(runs, key=lambda item: (item[0], str(item[1])))[1:]


def subplot_row(fig, default_ax, n):
    default_ax.remove()
    fig.set_size_inches(4.7 * n, 4.5)
    return fig.subplots(1, n, squeeze=False)[0]


def colorbar(fig, axes, image, matrices, limits, label):
    low, high = limits
    tails = (any(m.min() < low for m in matrices), any(m.max() > high for m in matrices))
    extend = {(False, False): "neither", (True, False): "min", (False, True): "max", (True, True): "both"}[tails]
    fig.colorbar(image, ax=axes, fraction=0.045, pad=0.035, label=label, extend=extend)


def image(ax, matrix, title, labels, cmap, limits):
    ax.set(title=title, xlabel=labels, ylabel=labels)
    return ax.imshow(matrix, origin="upper", aspect="equal", interpolation="nearest",
                     cmap=cmap, vmin=limits[0], vmax=limits[1])


def plot_sc(out, sc_model, model):
    # The MFM takes the subject's empirical SC as a fixed input. Both panels
    # show that same matrix in the actual units used for integration.
    normalization = model["sc_normalization"]
    limits = (0, float(sc_model.max()) if normalization == "as_supplied" else 1)
    label = {"column_max": "Poids / maximum de sa colonne",
             "global_max": "Poids / maximum global",
             "as_supplied": model["sc_model_weights"]}[normalization]

    def draw(fig, ax):
        axes = subplot_row(fig, ax, 2)
        fig.suptitle("SC du sujet utilisée comme entrée fixe du MFM")
        for panel, title in zip(axes, ("SC empirique · unités du modèle", "SC injectée dans le MFM")):
            heat = image(panel, sc_model, title, "Nœuds", "turbo", limits)
        colorbar(fig, axes, heat, (sc_model,), limits, label)
    io.figure(out, "01_sc_du_sujet", draw)
    return limits


def fisher_z(r):
    z = np.zeros_like(r)
    z[io.EDGE] = np.arctanh(np.clip(r[io.EDGE], -1 + 1e-12, 1 - 1e-12))
    return z + z.T


def plot_fc(out, observed, simulated, limits, empirical_label):
    difference = simulated - observed
    bound = max(float(np.max(np.abs(difference))), np.finfo(float).eps)

    def draw(fig, ax):
        axes = subplot_row(fig, ax, 2)
        fig.suptitle("FC empirique et simulée · échelle commune en z de Fisher")
        for panel, matrix, title in zip(axes, (observed, simulated), (empirical_label, "MFM simulé")):
            heat = image(panel, matrix, title, "Nœuds", "turbo", limits)
        colorbar(fig, axes, heat, (observed, simulated), limits, "z de Fisher")
    io.figure(out, "02_fc_statique", draw)
    def draw_difference(fig, ax):
        heat = image(ax, difference, "FC · MFM − empirique", "Nœuds", "RdBu_r", (-bound, bound))
        colorbar(fig, ax, heat, (difference,), (-bound, bound), "Écart de FC · Δz")
    io.figure(out, "02_fc_ecart", draw_difference)


def fcd_panel(out, observed, simulated, version):
    def draw(fig, ax):
        axes = subplot_row(fig, ax, 2)
        fig.suptitle("FCD · 10 volumes · pas de 1 volume")
        limits = (0, 0.5)
        for panel, matrix, label in zip(axes, (observed, simulated), (f"Empirique {version}", "MFM simulé")):
            heat = image(panel, matrix, label, "Fenêtre", "jet", limits)
        colorbar(fig, axes, heat, (observed, simulated), limits, "Corrélation entre FC")
    io.figure(out, "03_fcd_10_volumes", draw)


def density_panel(out, observed, simulated, version):
    def draw(fig, ax):
        axes = subplot_row(fig, ax, 2)
        fig.suptitle("Distribution de la récurrence et de la vitesse dFC")
        for i, (panel, label, xmin) in enumerate(zip(axes, ("χ · récurrence", "ν · vitesse"), (-1, 0))):
            samples = (observed[i], simulated[i])
            lower = max(xmin, min(values.min() for values in samples) - 0.08)
            upper = min(xmin + 2, max(values.max() for values in samples) + 0.08)
            support = np.linspace(lower, upper, 350)
            for values, name, color in zip(samples, (f"Empirique {version}", "MFM simulé"), (EMP_COLOR, MFM_COLOR)):
                description = f"{name} (n = {len(values)})"
                if len(values) > 1 and np.ptp(values) > 1e-12:
                    panel.plot(support, gaussian_kde(values)(support), color=color, lw=2, label=description)
                else:
                    panel.axvline(values[0], color=color, label=description)
            panel.set(xlabel=label, ylabel="Densité", xlim=(lower, upper), ylim=(0, None))
            panel.legend(frameon=False)
    io.figure(out, "04_distributions_dfc", draw)


def quantile_panel(out, observed, simulated, version):
    def draw(fig, ax):
        axes = subplot_row(fig, ax, 2)
        fig.suptitle("Cibles 2026 · quantiles empiriques et simulés")
        for panel, begin, label in zip(axes, (0, 5), ("χ · récurrence", "ν · vitesse")):
            section = slice(begin, begin + 5)
            panel.plot(PERCENTILES, observed[section], "o-", color=EMP_COLOR, label=f"Empirique {version}")
            panel.plot(PERCENTILES, simulated[section], "o-", color=MFM_COLOR, label="MFM simulé")
            panel.set(xlabel="Percentile", ylabel=label, xticks=PERCENTILES)
            panel.legend(frameon=False)
    io.figure(out, "05_cibles_2026", draw)


def static_comparison(out, cfg, model, r_emp, r_sim, sc_model):
    version = cfg["version"]
    fc_source = "BOLD retenu"
    if cfg.get("fc_reference_z") and cfg["censor_mode"] == "none":
        # Le panneau de référence reprend exactement l'aperçu fc_v4_adaptee.
        ref_path = io.source(cfg, "fc_reference_z")
        z_emp = io.array_file(ref_path).copy()
        if (z_emp.shape != sc_model.shape or not np.isfinite(z_emp[io.EDGE]).all()
                or not np.allclose(z_emp[io.EDGE], z_emp.T[io.EDGE], atol=1e-8)):
            raise ValueError(f"FC de référence incompatible : {ref_path}")
        np.fill_diagonal(z_emp, 0)
        fc_source = ref_path.name
        emp_label = f"Empirique {version} · FC fournie"
    else:
        z_emp = fisher_z(r_emp)
        emp_label = f"Empirique {version} · BOLD retenu"
    z_sim = fisher_z(r_sim)
    fc_limits = fc_display_limits(z_emp)
    sc_limits = plot_sc(out, sc_model, model)
    plot_fc(out, z_emp, z_sim, fc_limits, emp_label)
    scan = f"{cfg['subject']} / {cfg['session']}"
    for name, matrix, label in (("fc_empirique", z_emp, emp_label), ("fc_simulee", z_sim, "MFM simulé")):
        plot_static_matrix(out, name, matrix, f"{label} · échelle commune {version} · {scan}",
                           fc_limits, f"z de Fisher · bornes empiriques {version}")
    return fc_source, fc_limits, sc_limits


def render_comparison(out, cfg, model, x, y, grid, sc, sc_model):
    """Tracer les matrices sur la même grille de volumes retenus et renvoyer les cibles 2026."""
    if x.shape != y.shape or len(grid) != len(x):
        raise ValueError("Les séries empirique et simulée doivent partager la même grille.")
    data = out / "data"
    tr = cfg["tr_seconds"]
    emp, _, _, _, chi_emp, nu_emp, _ = feature_set(x, grid, tr)
    sim, _, _, fcd_targets, chi_sim, nu_sim, _ = feature_set(y, grid, tr)
    fcd_emp = fcd_ten_volumes(x, grid)
    fcd_sim = fcd_ten_volumes(y, grid)
    r_emp, r_sim = fc(x), fc(y)
    fc_source, fc_limits, sc_limits = static_comparison(out, cfg, model, r_emp, r_sim, sc_model)
    edges = r_emp[io.EDGE], r_sim[io.EDGE]
    diff = edges[1] - edges[0]
    pattern = float(np.corrcoef(*edges)[0, 1]) if all(np.ptp(v) > 0 for v in edges) else None
    width = volumes(60, tr)
    summary = {"volumes_retenus": len(grid), "tr_secondes": tr, "fenetre_cibles_volumes": width,
               "fenetre_fcd_figure_volumes": 10 if fcd_emp is not None else None,
               "fcd_figure_shape": list(fcd_emp.shape) if fcd_emp is not None else None,
               "fcd_cibles_shape": list(fcd_targets.shape), "n_vitesses_mfm": len(nu_sim),
               "fc_bornes_affichage_z_empirique": list(fc_limits),
               "fc_source_empirique_figure": fc_source, "fc_palette": "turbo",
               "fc_source_metriques_pearson": "BOLD empirique et simulé sur la même grille",
               "sc_source_empirique": model["sc_model"], "sc_bornes_affichage_empirique": list(sc_limits),
               "normalisation_sc": model["sc_normalization"], "fc_correlation_entre_liens": pattern,
               "fc_mae_pearson": float(np.mean(np.abs(diff))), "fc_rmse_pearson": float(np.sqrt(np.mean(diff**2))),
               "nu_mediane_empirique": float(np.median(nu_emp)), "nu_mediane_mfm": float(np.median(nu_sim)),
               "distance_euclidienne_cibles": float(np.linalg.norm(sim - emp))}
    io.write_json(data / "comparaison.json", summary)
    io.write_table(data / "comparaison_cibles.csv", ["cible", "empirique", "mfm", "mfm_moins_empirique"],
                   zip(TARGETS, emp, sim, sim - emp))
    np.save(data / "fcd_simulee_cibles_60s.npy", fcd_targets, allow_pickle=False)
    scan = f"{cfg['subject']} / {cfg['session']}"
    version = cfg["version"]
    if fcd_emp is not None:
        for name, matrix, label in (("fcd_empirique", fcd_emp, f"empirique {version}"),
                                    ("fcd_simulee", fcd_sim, "simulée MFM")):
            np.save(data / f"{name}_10vol.npy", matrix, allow_pickle=False)
            io.heatmap(out, name, matrix, f"FCD {label} · {scan} · 10 volumes", "Fenêtre",
                       style="fcd_article", limits=(0, 0.5))
        fcd_panel(out, fcd_emp, fcd_sim, version)
    density_panel(out, (chi_emp, nu_emp), (chi_sim, nu_sim), version)
    quantile_panel(out, emp, sim, version)
    return emp, sim


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", help="Identifiant ou dossier de la simulation. Par défaut, la dernière simulation terminée du pilote.")
    parser.add_argument("--run", help="Identifiant du lancement de sortie")
    parser.add_argument("--retracer-sc-fc", metavar="RUN", help="Remplacer les figures statiques dans une comparaison terminée")
    args = parser.parse_args()
    if args.retracer_sc_fc:
        if args.source_run or args.run:
            parser.error("--retracer-sc-fc utilise la source enregistrée dans la comparaison.")
        tag = "_".join([io.token(args.retracer_sc_fc), *PILOTE])
        out = io.stage_dir("10_jumeau_comparaison") / tag
        comparison = io.read_json(out / "manifest.json")
        if comparison["status"] != "complete":
            parser.error("La comparaison doit être terminée avant de remplacer ses figures.")
        args.source_run = str(Path(comparison["inputs"]["manifest_controle"]["path"]).parent)
    manifest_path, state, cfg, model, seed, bold_simulated = control_manifest(args.source_run)
    x, grid, n, sc, names, files = io.load_inputs(cfg)
    sc_model, model_files = setup(cfg, model, names)
    files.update(model_files)
    for name, path in files.items():
        if io.digest(path) != state["inputs"][name]["sha256"]:
            raise ValueError(f"Entrée modifiée depuis la simulation : {name}")
    if cfg.get("fc_reference_z") and cfg["censor_mode"] == "none":
        path = io.source(cfg, "fc_reference_z")
        old_ref = state["inputs"].get("fc_reference_z")
        if old_ref and io.digest(path) != old_ref["sha256"]:
            raise ValueError("La FC de référence a été modifiée depuis la simulation.")
        files["fc_reference_z"] = path
    y = np.load(bold_simulated, allow_pickle=False)
    if y.shape != (n, len(names)) or not np.isfinite(y).all():
        raise ValueError("Le BOLD simulé et la grille empirique sont incompatibles.")
    if args.retracer_sc_fc:
        for key, path in dict(files, bold_simule=bold_simulated, manifest_controle=manifest_path).items():
            if io.digest(path) != comparison["inputs"][key]["sha256"]:
                raise ValueError(f"Entrée différente de la comparaison enregistrée : {key}")
        dest = out / "figures"
        with TemporaryDirectory(dir=dest) as tmp:
            staging = Path(tmp)
            (staging / "figures").mkdir()
            _, _, sc_limits = static_comparison(staging, cfg, model, fc(x), fc(y[grid]), sc_model)
            summary_path = out / "data/comparaison.json"
            summary = io.read_json(summary_path)
            summary.update(sc_source_empirique=model["sc_model"], sc_bornes_affichage_empirique=list(sc_limits))
            io.write_json(staging / "comparaison.json", summary)
            for path in sorted((staging / "figures").iterdir()):
                path.replace(dest / path.name)
                print(dest / path.name)
            (staging / "comparaison.json").replace(summary_path)
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = io.run_dir(cfg, "10_jumeau_comparaison", args.run or stamp)
    settings = {"run_controle": manifest_path.parent.name, "pilote": cfg,
                "parametres_mfm": model["control_parameters"], "graine": seed}
    files.update(bold_simule=bold_simulated, manifest_controle=manifest_path, script_comparaison=__file__)
    with io.record(out, "Comparaison empirique V4 et MFM simulé", settings, files):
        render_comparison(out, cfg, model, x, y[grid], grid, sc, sc_model)
    print(f"Source : {manifest_path.parent}")
    print(out / "figures/02_fc_statique.png")
    print(out / "figures/03_fcd_10_volumes.png")


if __name__ == "__main__":
    main()
