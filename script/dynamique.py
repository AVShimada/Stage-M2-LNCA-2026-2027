"""FC, dynamique par fenêtres, méta-connectivité et statistiques de validation."""
import math
import numpy as np

from commun import EDGE, SPEED_STYLES, figure, heatmap, write_json, write_table

PERCENTILES = (5, 25, 50, 75, 95)
TARGETS = [f"{v}_p{p}" for v in ("chi", "nu") for p in PERCENTILES]


def unit(x, axis):
    z = x - x.mean(axis=axis, keepdims=True)
    norm = np.linalg.norm(z, axis=axis, keepdims=True)
    if np.any(norm == 0) or not np.isfinite(norm).all():
        raise ValueError("Corrélation indéfinie : série, région ou vecteur de liens constant")
    return z / norm


def fc(x):
    z = unit(x, 0)
    return np.clip(z.T @ z, -1, 1)


def blocks(grid):
    return np.split(np.arange(len(grid)), np.flatnonzero(np.diff(grid) != 1) + 1)


def volumes(seconds, tr):
    return max(1, int(math.floor(seconds / tr + 0.5)))


def windows(x, grid, width, step):
    if width < 3:
        raise ValueError("Une fenêtre doit contenir au moins trois volumes")
    starts = [i for block in blocks(grid) if len(block) >= width
              for i in range(block[0], block[-1] - width + 2, step)]
    if len(starts) < 2:
        raise ValueError(f"Fenêtres continues insuffisantes : largeur {width} TR, pas {step} TR")
    return np.array([fc(x[i:i+width])[EDGE] for i in starts]), grid[starts]


def paired_speeds(stream, starts, width, step, rng=None):
    # Within one contiguous segment, compare frames separated by one window width.
    first = np.arange(len(starts) - width // step)
    second = first + width // step
    ok = starts[second] - starts[first] == width
    first, second = first[ok], second[ok]
    if not len(first):
        raise ValueError(f"Aucune paire de fenêtres contiguës de {width} TR")
    z = unit(stream, 1)
    if rng is not None:
        z = z[rng.permutation(len(z))]
    return 1 - np.clip(np.einsum("ij,ij->i", z[first], z[second]), -1, 1), starts[first]


def feature_set(x, grid, tr):
    width, step = volumes(60, tr), volumes(2, tr)
    stream, starts = windows(x, grid, width, step)
    z = unit(stream, 1)
    fcd = np.clip(z @ z.T, -1, 1)
    chi = fcd[np.triu_indices(len(fcd), 1)]
    speeds, details = [], []
    for w in range(math.ceil(55 / tr), math.floor(65 / tr) + 1):
        s, t = windows(x, grid, w, w)
        nu, t0 = paired_speeds(s, t, w, w)
        speeds.extend(nu)
        details.extend(zip([w] * len(nu), t0, nu))
    if not speeds:
        raise ValueError("La grille 55–65 s ne contient aucune largeur réalisable au TR indiqué")
    nu = np.array(speeds)
    target = np.r_[np.percentile(chi, PERCENTILES), np.percentile(nu, PERCENTILES)]
    return target, stream, starts, fcd, chi, nu, details


def mc_matrix(stream):
    z = unit(stream, 0).astype("float32")
    return np.clip(z.T @ z, -1, 1)


def validation(x, sc, stream, mc):
    scvec = unit(sc[EDGE], 0)
    lam = np.clip(unit(stream, 1) @ scvec, -1, 1)
    # 4,950 unique links; sample the unique pairs of links exactly once.
    values = {"rho": fc(x)[EDGE], "lambda": lam, "MC": mc[np.triu_indices(len(mc), 1)]}
    return {f"{k}_p{p}": float(q) for k, v in values.items()
            for p, q in zip((5, 50, 95), np.percentile(v, [5, 50, 95]))}, lam


def static(out, x, sc, names):
    r = fc(x)
    np.save(out / "data/fc_r.npy", r)
    np.save(out / "data/sc.npy", sc)
    if np.any(np.abs(r[EDGE]) >= 1):
        raise ValueError("FC parfaite hors diagonale : vérifier les séries avant la transformation de Fisher")
    z = np.zeros_like(r)
    z[EDGE] = np.arctanh(r[EDGE])
    z += z.T
    np.save(out / "data/fc_fisher_z.npy", z)
    metrics = {}
    for key, m in (("fc_r", r), ("fc_z", z), ("sc", sc)):
        v = m[EDGE]
        weights = m.sum(axis=1) - np.diag(m)
        metrics.update({f"{key}_{k}": float(vv) for k, vv in
                        (("mean", v.mean()), ("std", v.std()), ("median", np.median(v)), ("sum", v.sum()),
                         ("abs_mean", np.abs(v).mean()), ("nonzero_fraction", (v != 0).mean()),
                         ("mean_strength", weights.mean()), ("std_strength", weights.std()))})
    write_json(out / "data/metrics.json", metrics)
    write_table(out / "data/edges.csv", ["edge_0", "roi_i", "roi_j"],
                ((k, names[i], names[j]) for k, (i, j) in enumerate(zip(*EDGE))))
    heatmap(out, "fc", r, "FC statique · corrélation de Pearson")
    heatmap(out, "sc", sc, "Connectivité structurelle · poids fournis", signed=False)


def phase_random(x, grid, rng):
    result = x.copy()
    for block in blocks(grid):
        n = len(block)
        ft = np.fft.rfft(x[block], axis=0)
        phase = rng.uniform(-np.pi, np.pi, len(ft))
        phase[0] = 0
        if n % 2 == 0:
            phase[-1] = 0
        result[block] = np.fft.irfft(ft * np.exp(1j * phase[:, None]), n=n, axis=0)
    return result


def m1(out, x, grid, tr, cfg):
    data = out / "data"
    w = volumes(20, tr)
    stream, starts = windows(x, grid, w, 1)
    z = unit(stream, 1)
    fcd = np.clip(z @ z.T, -1, 1)
    np.save(data / "stream_20s.npy", stream)
    np.save(data / "starts_20s.npy", starts)
    np.save(data / "fcd_20s.npy", fcd)
    np.save(data / "mc_20s.npy", mc_matrix(stream))
    speed, t = paired_speeds(stream, starts, w, 1)
    write_table(data / "speeds_20s.csv", ["first_volume_0", "speed"], zip(t, speed))
    heatmap(out, "fcd_20s", fcd, "FCD · fenêtres de 20 s, pas de 1 TR", "Fenêtre retenue", style="fcd_m1")
    mc_view = np.load(data / "mc_20s.npy", mmap_mode="r")[::5, ::5]
    heatmap(out, "mc_20s", mc_view, "Méta-connectivité · un lien affiché sur cinq", "Indice de lien (× 5)",
            style="mc_m1", limits=(0, max(np.finfo(float).eps, float(np.percentile(mc_view, 99)))))
    settings = cfg["m1"]
    strategy = settings["speed_definition"]
    if strategy not in ("pooled_sliding", "median_nonoverlap"):
        raise ValueError("m1.speed_definition inconnue")
    limit = settings["max_volumes_for_speed"]
    if limit is not None:
        keep = grid < limit
        x, grid = x[keep], grid[keep]
    rng = np.random.default_rng(cfg["seed"])
    rows, stats, skipped = [], {"Vtyp_20s": float(np.median(speed))}, []
    ranges = {"mid": (15, 60), "long": (60, 210)}
    n_sur = settings["surrogates"]
    if type(n_sur) is not int or n_sur < 0:
        raise ValueError("m1.surrogates doit être un entier positif ou nul")
    for rep in range(n_sur + 1):
        signal = x if rep == 0 else phase_random(x, grid, rng)
        for name, (low, high) in ranges.items():
            collected = {}
            for width in range(round(low / tr), round(high / tr) + 1, 10):
                if width < 3 or not any(len(b) >= 2 * width for b in blocks(grid)):
                    skipped.append([name, width])
                    continue
                step = 1 if strategy == "pooled_sliding" else width
                s, start = windows(signal, grid, width, step)
                variants = [("emp", s, None)] if rep == 0 else [("phase", s, None)]
                if rep:
                    se, te = windows(x, grid, width, step)
                    if not np.array_equal(te, start):
                        raise ValueError("Fenêtres incompatibles")
                    variants.append(("shuffle", se, rng))
                for kind, sample, gen in variants:
                    v, ti = paired_speeds(sample, start, width, step, gen)
                    rows.extend((name, rep, kind, width, int(i), float(val)) for i, val in zip(ti, v))
                    collected.setdefault(kind, []).extend(v if strategy == "pooled_sliding" else [np.median(v)])
            if not collected:
                raise ValueError(f"Aucune taille de fenêtre exploitable pour Vtyp_{name}")
            for kind, v in collected.items():
                stats[f"Vtyp_{name}_{kind}_{rep}"] = float(np.median(v))
    write_table(data / "speeds_multiscale.csv", ["range","rep","type","width_TR","first_volume_0","nu"], rows)
    write_json(data / "metrics.json", stats)
    write_json(data / "omissions.json", sorted(set(tuple(v) for v in skipped)))
    def draw(fig, ax):
        for kind in ("emp", "phase", "shuffle"):
            v = [r[-1] for r in rows if r[0] == "mid" and r[2] == kind]
            if v:
                color, line = SPEED_STYLES[kind]
                ax.hist(v, bins=np.linspace(0,2,41), density=True, histtype="step", label=kind,
                        color=color, linestyle=line, linewidth=2)
        ax.set(xlabel="Vitesse (1 − r)", ylabel="Densité", title="Vitesses dFC · fenêtres moyennes M1")
        ax.legend()
    figure(out, "vitesses", draw, dark=True)


def empirical_targets(out, x, grid, tr, sc):
    data = out / "data"
    target, stream, starts, fcd, chi, nu, details = feature_set(x, grid, tr)
    mc = mc_matrix(stream)
    checks, lam = validation(x, sc, stream, mc)
    for name, mat in (("stream", stream), ("fcd", fcd), ("starts_TR", starts), ("mc", mc)):
        np.save(data / f"{name}.npy", mat)
    write_table(data / "speeds.csv", ["width_TR", "first_volume_0", "nu"], details)
    write_table(data / "coupling.csv", ["first_volume_0", "lambda"], zip(starts, lam))
    write_json(data / "targets.json", dict(zip(TARGETS, target.tolist())))
    write_json(data / "validation.json", checks)
    write_json(data / "windows.json", dict(width_TR=volumes(60,tr), step_TR=volumes(2,tr),
               width_s=tr*volumes(60,tr), step_s=tr*volumes(2,tr), n_chi=len(chi), n_nu=len(nu)))
    heatmap(out, "fcd", fcd, "FCD · profil 2026", "Fenêtre retenue", style="fcd_article")
    def draw(fig, ax):
        for vals, label, line in ((chi, "χ · récurrence", "-"), (nu, "ν · vitesse", "--")):
            ax.hist(vals, bins=40, density=True, histtype="step", label=label, color="black", linestyle=line)
        ax.set(xlabel="Valeur", ylabel="Densité", title="Distributions utilisées pour ajuster le jumeau")
        ax.legend()
    figure(out, "distributions", draw)

