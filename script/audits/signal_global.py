"""Signal moyen des 100 régions standardisées et FC avant/après sa régression."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import commun as io


def regional_series(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or x.shape[1] != 100 or x.shape[0] < 3 or not np.isfinite(x).all():
        raise ValueError(f"Série attendue : temps × 100 régions, valeurs finies, reçu {x.shape}.")
    x = x - x.mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    if np.any(sd <= 0):
        raise ValueError("Région constante dans la série.")
    z = x / sd
    return z.T @ z / (len(z) - 1), z.mean(axis=1)


def project_global(matrix):
    r = np.asarray(matrix, dtype=float)
    if (r.shape != (100, 100) or not np.isfinite(r).all()
            or not np.allclose(r, r.T, atol=1e-7) or not np.allclose(np.diag(r), 1.0)
            or np.linalg.eigvalsh(r)[0] < -1e-5):
        raise ValueError("Matrice de corrélations de Pearson invalide.")
    # c_i = cov(z_i, g), g est la moyenne régionale de z.
    c = r.mean(axis=1)
    var_g = float(c.mean())
    if var_g < -1e-12:
        raise ValueError("Variance du signal moyen négative.")
    if var_g <= 1e-12:
        return r.copy(), dict(variance_signal_moyen=0., r2_regional_moyen=0.)
    residual = r - np.outer(c, c) / var_g
    if np.any(np.diag(residual) <= 1e-12):
        raise ValueError("Variance régionale résiduelle nulle après régression.")
    sd = np.sqrt(np.diag(residual))
    return residual / np.outer(sd, sd), dict(variance_signal_moyen=var_g,
                                            r2_regional_moyen=float(np.mean(c**2) / var_g))


def link_metrics(a, b):
    x, y = a[io.EDGE], b[io.EDGE]
    return dict(rmse=float(np.sqrt(np.mean((y - x)**2))),
                biais_moyen=float(np.mean(y - x)),
                correlation_profils=float(np.corrcoef(x, y)[0, 1]) if min(np.ptp(x), np.ptp(y)) > 0 else None)


def draw_summary(out, stats, comparisons, subject, session, seed):
    def draw(fig, ax):
        ax.remove()
        fig.set_size_inches(10, 4.9)
        left, right = fig.subplots(1, 2)
        values = [100 * row["r2_regional_moyen"] for row in stats.values()]
        bars = left.bar(range(len(values)), values, color="#277481", width=.65)
        left.set_xticks(range(len(values)), [s.replace(" ", "\n", 1) for s in stats])
        left.set(title="Variance régionale expliquée", ylabel="R² moyen (%)", ylim=(0, max(5, max(values) * 1.25)))
        left.bar_label(bars, fmt="%.1f", padding=4)
        pos = np.arange(len(comparisons))
        for offset, key, title, color in [(-.18, "avant", "Avant", "#57697b"), (.18, "apres", "Après", "#cf7447")]:
            vals = [v[key]["rmse"] for v in comparisons.values()]
            bars = right.bar(pos + offset, vals, .35, label=title, color=color)
            right.bar_label(bars, fmt="%.3f", padding=3, fontsize=9)
        right.set_xticks(pos, [f"{a}\n↔ {b}" for a, b in comparisons])
        right.set(title="Écart de FC avant et après régression", ylabel="RMSE des 4 950 liens (r)")
        right.margins(y=.2)
        right.legend(frameon=False)
        for axis in (left, right):
            axis.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"Signal moyen des 100 régions · {subject} / {session} · jumeau, graine {seed}")
    io.figure(out, "signal_global_fc", draw)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilote", required=True, help="Identifiant UTC du run pilote V1/V4")
    p.add_argument("--jumeau", required=True, help="Identifiant UTC du run d'ajustement conjoint FC/dFC")
    p.add_argument("--racine-output", type=Path, default=io.OUTPUT)
    p.add_argument("--depuis-matrices", action="store_true", help="Analyser les exports téléchargés")
    args = p.parse_args()
    root = args.racine_output.resolve()
    subject, session, version = "sub-0001", "ses-1", "V4"
    suffix = f"_{subject}_{session}_{version}"
    pilot = io.token(args.pilote) + suffix
    twin = root / io.STAGE_PATHS["10_ajuster_mfm"] / "fc_dfc" / (io.token(args.jumeau) + suffix)
    pilot_dir = root / io.STAGE_PATHS["02_comparaison_fcd_v1_v4"] / pilot
    files = dict(manifeste_pilote=pilot_dir / "manifest.json", manifeste_jumeau=twin / "manifest.json",
                 resultat_jumeau=twin / "data/resultat.json", simulation=twin / "data/bold_simule_volumes_retenus.npy",
                 script=Path(__file__).resolve())
    pm, jm = (io.read_json(files[k]) for k in ("manifeste_pilote", "manifeste_jumeau"))
    cfg = pm["configuration"]
    if (pm["status"] != "complete" or jm["status"] != "complete"
            or pm["inputs"]["V4"]["sha256"] != jm["inputs"]["bold"]["sha256"]
            or cfg["censor_mode"] != "none"):
        raise ValueError("Les deux runs doivent être terminés et appariés sur le même BOLD V4 consécutif.")
    result = io.read_json(files["resultat_jumeau"])
    seed = result["graine_figures"]
    sim = io.array_file(files["simulation"])
    if sim.shape != (pm["n_volumes_per_version"], 100):
        raise ValueError(f"Grilles pilote/jumeau incompatibles : {sim.shape}")
    sim_fc, sim_global = regional_series(sim)
    raw, traces = {}, {}
    if args.depuis_matrices:
        preview = root / io.STAGE_PATHS["01_apercu_sc_fc"] / pilot
        static = root / io.STAGE_PATHS["01_sc_fc"] / pilot
        files.update(apercu=preview / "manifest.json", statique=static / "manifest.json")
        a, s = io.read_json(files["apercu"]), io.read_json(files["statique"])
        if (a["status"] != "complete" or s["status"] != "complete"
                or s["inputs"]["bold"]["sha256"] != pm["inputs"]["V4"]["sha256"]
                or any(a["configuration"][k] != cfg[k] for k in ("subject", "session", "atlas", "version"))):
            raise ValueError("Aperçu ou FC statique incompatible avec le pilote.")
        for v in ("V1", "V4"):
            label = f"{v} fournie"
            files[label] = preview / "data" / f"fc_{v.lower()}.npy"
            r = np.tanh(io.array_file(files[label]))
            np.fill_diagonal(r, 1.)
            raw[label] = r
        files["V4 BOLD"] = static / "data/fc_r.npy"
        raw["V4 BOLD"] = io.array_file(files["V4 BOLD"])
        pairs = [("V1 fournie", "V4 fournie"), ("V4 BOLD", "Jumeau")]
    else:
        for v in ("V1", "V4"):
            path = Path(pm["inputs"][v]["path"])
            if io.digest(path) != pm["inputs"][v]["sha256"]:
                raise ValueError(f"Série {v} différente de celle du pilote.")
            files[v] = path
            x = io.array_file(path)
            if cfg["bold_orientation"] == "regions_time":
                x = x.T
            elif cfg["bold_orientation"] != "time_regions":
                raise ValueError("Orientation BOLD inconnue.")
            if x.shape != sim.shape:
                raise ValueError(f"Dimensions BOLD {v} inattendues : {x.shape}")
            raw[v], traces[v] = regional_series(x)
        pairs = [("V1", "V4"), ("V4", "Jumeau")]
    raw["Jumeau"] = sim_fc
    removed, stats = {}, {}
    for label, matrix in raw.items():
        removed[label], stats[label] = project_global(matrix)
        stats[label].update(mediane_fc_avant=float(np.median(matrix[io.EDGE])),
                            mediane_fc_apres=float(np.median(removed[label][io.EDGE])))
    comparisons = {(a,b): dict(avant=link_metrics(raw[a], raw[b]), apres=link_metrics(removed[a], removed[b]))
                   for a,b in pairs}
    mode = "matrices" if args.depuis_matrices else "series"
    settings = dict(sujet=subject, session=session, mode=mode, graine=seed, run_pilote=pilot,
                    run_jumeau=twin.name, parametres={k:result["retenu"][k] for k in ("G","w","sigma")},
                    signal_moyen="Moyenne non pondérée des 100 régions centrées et réduites dans chaque série.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = root / io.STAGE_PATHS["12_evaluer_jumeaux"] / "signal_global" / (stamp + suffix)
    with io.record(out, "Signal global et écarts de FC", settings, files) as data:
        summary = dict(mode=mode, n_volumes=len(sim), n_liens=len(io.EDGE[0]), graine_jumeau=seed,
                       series=stats, comparaisons={f"{a} / {b}": row for (a,b),row in comparisons.items()})
        if traces:
            a, b = traces["V1"], traces["V4"]
            summary["correlation_signaux_moyens_V1_V4"] = (float(np.corrcoef(a,b)[0,1])
                    if min(a.var(), b.var()) > 1e-12 else None)
            io.write_table(data / "signaux_moyens.csv", ["volume", "temps_s", "V1", "V4", "jumeau"],
                           zip(range(len(sim)), np.arange(len(sim))*cfg["tr_seconds"], a, b, sim_global))
        draw_summary(out, stats, comparisons, subject, session, seed)
        io.write_json(data / "resume.json", summary)
    print(out / "data/resume.json")


if __name__ == "__main__":
    main()
