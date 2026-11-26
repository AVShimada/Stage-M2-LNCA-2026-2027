"""Tracer les figures du pilote et retracer un lancement déjà terminé."""

import argparse
import csv
from datetime import datetime, timezone

import numpy as np
from scipy.stats import gaussian_kde

from commun import EDGE, stage_dir, SPEED_STYLES, figure, heatmap, load_inputs, read_json


FC_LIMITS_M1 = (-0.40065664835678494, 1.2926806413667378)
SC_LIMITS_M1 = (-0.08296555948553058, 4.607515607717042)
SPEED_LABELS = {"emp": "Empirique", "phase": "Phase randomisée", "shuffle": "Fenêtres mélangées"}


def fc_display_limits(z):
    low, high = np.percentile(z[EDGE], [2, 98])
    return (float(low), float(high)) if low < high else FC_LIMITS_M1


def sc_display_limits(sc, weights):
    weights = weights.casefold()
    if "sift2" in weights and "cm_log" in weights:
        return SC_LIMITS_M1
    low, high = float(sc.min()), float(sc.max())
    return (low, high) if low < high else (min(0.0, low), max(1.0, high))


def plot_static_matrix(out, name, matrix, title, limits, label):
    display = matrix.copy()
    np.fill_diagonal(display, 0)
    low, high = limits
    below = display[EDGE].min() < low
    above = display[EDGE].max() > high
    extend = {(False, False): "neither", (True, False): "min", (False, True): "max", (True, True): "both"}

    def draw(fig, ax):
        im = ax.imshow(display, cmap="turbo", origin="upper", aspect="equal", interpolation="nearest",
                       vmin=low, vmax=high)
        ax.set(title=title, xlabel="Nœuds", ylabel="Nœuds")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend=extend[below, above], label=label)

    figure(out, name, draw)


def plot_static_fc(out, z):
    plot_static_matrix(out, "fc", z, "FC statique · z de Fisher · échelle adaptée",
                       fc_display_limits(z), "z de Fisher")


def plot_static_sc(out, sc, weights):
    plot_static_matrix(out, "sc", sc, f"SC · {weights}", sc_display_limits(sc, weights), weights)


def density_figure(out, name, series, limits, xlabel, title, *, dark=False):
    x = np.linspace(*limits, 400)

    def draw(fig, ax):
        for values, label, color, linestyle in series:
            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]
            if not len(values):
                continue
            options = dict(color=color, linestyle=linestyle, linewidth=2, label=f"{label} (n = {len(values):,})".replace(",", " "))
            if len(values) > 1 and np.ptp(values) > 0:
                try:
                    density = gaussian_kde(values)(x)
                except np.linalg.LinAlgError:
                    pass
                else:
                    ax.plot(x, density, **options)
                    continue
            ax.axvline(values[0], **options)
        ax.set(xlim=limits, ylim=(0, None), xlabel=xlabel, ylabel="Densité", title=title)
        if ax.lines:
            ax.legend()

    figure(out, name, draw, dark=dark)


def plot_m1_speeds(out, speeds):
    samples = {kind: [] for kind in SPEED_LABELS}
    for scale, kind, value in speeds:
        if scale == "mid":
            samples[kind].append(float(value))
    upper = max(1.2, *(max(values, default=0.0) for values in samples.values()))
    series = [(values, SPEED_LABELS[kind], *SPEED_STYLES[kind]) for kind, values in samples.items()]
    density_figure(out, "vitesses", series, (0, upper),
                   "Vitesse (1 − r)", "Vitesses dFC · fenêtres moyennes M1", dark=True)


def fcd_ten_volumes(x, grid):
    from dynamique import blocks, unit, windows

    if sum(max(0, len(block) - 9) for block in blocks(grid)) < 2:
        return None
    stream, _ = windows(x, grid, 10, 1)
    normalized = unit(stream, 1)
    return np.clip(normalized @ normalized.T, -1, 1)


def plot_targets(out, fcd_ten, chi, nu):
    if fcd_ten is not None:
        heatmap(out, "fcd", fcd_ten, "FCD · 10 volumes · pas de 1 volume", "Fenêtre",
                style="fcd_article", limits=(0, 0.5))
    series = [(chi, "χ · récurrence", "#242424", "-"), (nu, "ν · vitesse", "#0072BD", "--")]
    limits = (min(0.0, float(np.min(chi))), max(1.0, float(np.max(nu))))
    density_figure(out, "distributions", series, limits,
                   "Valeur", "Distributions empiriques · ajustement du jumeau")


def archive_old_fcd_figures(out, fcd_ten):
    names = ["fcd_cibles"] if fcd_ten is not None else ["fcd_cibles", "fcd"]
    old = [out / "figures" / f"{name}.{ext}" for name in names for ext in ("png", "pdf")]
    archive = out / "archives" / "figures_fcd" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    for path in old:
        if path.is_file():
            archive.mkdir(parents=True, exist_ok=True)
            path.rename(archive / path.name)


def retrace(run):
    from dynamique import volumes

    stages = ("01_sc_fc", "02_dynamique_m1", "03_cibles_2026")
    candidates = []
    for target in stage_dir(stages[-1]).glob("*/manifest.json"):
        tag = target.parent.name
        if run and tag != run and not tag.startswith(run + "_"):
            continue
        folders = [stage_dir(stage) / tag for stage in stages]
        manifests = [folder / "manifest.json" for folder in folders]
        if not all(path.is_file() for path in manifests):
            continue
        states = [read_json(path) for path in manifests]
        cfg = states[-1]["configuration"]
        if all(state["status"] == "complete" and state["configuration"] == cfg for state in states):
            candidates.append((states[-1]["started"], tag, cfg, folders))
    if not candidates:
        raise ValueError("Aucun pilote complet trouvé pour les étapes SC/FC, dynamique M1 et cibles 2026.")
    _, tag, cfg, (static, m1, targets) = max(candidates)

    z = np.load(static / "data/fc_fisher_z.npy", allow_pickle=False)
    sc = np.load(static / "data/sc.npy", allow_pickle=False)
    with (m1 / "data/speeds_multiscale.csv").open(newline="", encoding="utf-8-sig") as f:
        speeds = [(row["range"], row["type"], row["nu"]) for row in csv.DictReader(f)]
    fcd = np.load(targets / "data/fcd.npy", allow_pickle=False)
    visual_data = targets / "data/fcd_10_volumes.npy"
    if visual_data.is_file():
        fcd_ten = np.load(visual_data, allow_pickle=False)
    elif volumes(20, cfg["tr_seconds"]) == 10:
        fcd_ten = np.load(m1 / "data/fcd_20s.npy", allow_pickle=False)
    else:
        x, grid, *_ = load_inputs(cfg)
        fcd_ten = fcd_ten_volumes(x, grid)
    chi = fcd[np.triu_indices(len(fcd), 1)]
    with (targets / "data/speeds.csv").open(newline="", encoding="utf-8-sig") as f:
        nu = [float(row["nu"]) for row in csv.DictReader(f)]

    plot_static_fc(static, z)
    plot_static_sc(static, sc, cfg["sc_weights"])
    plot_m1_speeds(m1, speeds)
    plot_targets(targets, fcd_ten, chi, nu)
    if fcd_ten is not None and not visual_data.is_file():
        np.save(visual_data, fcd_ten, allow_pickle=False)
    archive_old_fcd_figures(targets, fcd_ten)
    print(f"Figures retracées pour {tag}")
    for folder in (static, m1, targets):
        print(folder / "figures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", help="Identifiant du lancement, dernier pilote complet par défaut")
    retrace(parser.parse_args().run)
