"""Adaptation individuelle des analyses de graphes réalisées en M1."""
import csv
import bct
import networkx as nx
import numpy as np

from commun import EDGE, figure, heatmap, write_json, write_table


def partition(weights, gamma, seed):
    graph = nx.from_numpy_array(weights)
    groups = (nx.community.louvain_communities(graph, weight="weight", resolution=gamma, seed=seed)
              if graph.number_of_edges() else [{i} for i in graph])
    ids = np.zeros(len(weights), dtype=int)
    for k, members in enumerate(groups, 1):
        ids[list(members)] = k
    return ids


def graph_metrics(out, names, outputs, settings, seed):
    metrics = {}
    for label, matrix in outputs.items():
        weights = np.maximum(matrix, 0)
        np.fill_diagonal(weights, 0)
        if weights.max() <= 0:
            raise ValueError(f"Graphe vide : {label}")
        parts = np.array([partition(weights, settings["resolution"], seed+i) for i in range(settings["repeats"])])
        agreement = np.mean(parts[:, :, None] == parts[:, None, :], axis=0)
        threshold = np.where(agreement >= settings["agreement_threshold"], agreement, 0)
        np.fill_diagonal(threshold, 0)
        ids = partition(threshold, 1.0, seed)
        strength = weights.sum(axis=1)
        frac = np.zeros(len(weights))
        for group in np.unique(ids):
            inside = weights[:, ids == group].sum(axis=1)
            frac += np.divide(inside, strength, out=np.zeros_like(inside), where=strength>0)**2
        participation = np.where(strength > 0, 1 - frac, 0)
        z = (strength - strength.mean()) / strength.std() if strength.std() else np.zeros_like(strength)
        np.savez(out / f"data/{label}_consensus.npz", partitions=parts, agreement=agreement, labels=ids)
        write_table(out / f"data/{label}_regions.csv", ["roi", "module", "strength", "strength_z", "participation"],
                    zip(names, ids, strength, z, participation))
        sort = np.argsort(ids, kind="stable")
        heatmap(out, label+"_modules", weights[np.ix_(sort, sort)], label+" · tri par consensus", signed=False)
        graph = nx.from_numpy_array(weights / weights.max())
        for i, j, attr in graph.edges(data=True):
            attr["distance"] = 1 / attr["weight"]
        total = sum(1/d for _, dist in nx.all_pairs_dijkstra_path_length(graph, weight="distance")
                    for d in dist.values() if d > 0)
        metrics[label+"_efficiency"] = total / (len(names)*(len(names)-1))
        metrics[label+"_clustering"] = nx.average_clustering(graph, weight="weight")
        def draw(fig, ax):
            ax.scatter(strength, participation, s=14, color="#30b7d4")
            ax.set(xlabel="Force pondérée", ylabel="Coefficient de participation", title=label+" · hubs")
            ax.xaxis.label.set_color("#30b7d4")
            ax.yaxis.label.set_color("#f27a1a")
        figure(out, label+"_hubs", draw, dark=True)
    write_json(out / "data/metrics.json", metrics)


def rich(out, fc, count, seed):
    if count < 1:
        raise ValueError("Il faut au moins un graphe nul")
    weights = np.maximum(fc, 0)
    np.fill_diagonal(weights, 0)
    obs = bct.rich_club_wu(weights)
    refs = []
    for i in range(count):
        null, _ = bct.null_model_und_sign(weights, bin_swaps=10, seed=seed+i)
        refs.append(bct.rich_club_wu(null))
    refs = np.array(refs)
    avg = refs.mean(axis=0)
    norm = np.divide(obs, avg, out=np.full_like(avg, np.nan), where=avg>0)
    k = np.arange(1, len(obs)+1)
    np.save(out / "data/rich_club_nulls.npy", refs)
    write_table(out / "data/rich_club.csv", ["k","phi","phi_null_mean","phi_normalized"], zip(k,obs,avg,norm))
    def draw(fig, ax):
        ax.plot(k, norm, color="white")
        ax.axhline(1, color="gray", linestyle="--")
        ax.set(xlabel="Seuil de degré k", ylabel="Coefficient normalisé", title="Rich club pondéré · FC positive")
    figure(out, "rich_club", draw, dark=True)


def distance(out, path, names, fc):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    key = "roi" if rows and "roi" in rows[0] else "label"
    coords = {r[key]: [float(r[k]) for k in ("x","y","z")] for r in rows}
    if len(rows) != 100 or set(coords) != set(names):
        raise ValueError("Le CSV des coordonnées doit associer roi,x,y,z aux 100 labels")
    pts = np.array([coords[n] for n in names])
    if not np.isfinite(pts).all():
        raise ValueError("Coordonnées non finies")
    dist = np.linalg.norm(pts[EDGE[0]] - pts[EDGE[1]], axis=1)
    v = fc[EDGE]
    keep = v >= 0
    dist, v = dist[keep], v[keep]
    bins = np.linspace(dist.min(), dist.max(), 9)
    which = np.searchsorted(bins[1:-1], dist, side="right")
    data = []
    for i in range(8):
        sel = which == i
        data.append([bins[i], bins[i+1], int(sel.sum()), v[sel].mean() if sel.any() else np.nan])
    write_table(out / "data/fc_distance.csv", ["from_mm","to_mm","n_edges","mean"], data)
    def draw(fig, ax):
        ax.plot((bins[:-1]+bins[1:])/2, [r[-1] for r in data], ".-", color="white")
        ax.set(xlabel="Distance euclidienne (mm)", ylabel="FC moyenne (z de Fisher)", title="Distance entre régions et FC positive")
    figure(out, "distance", draw, dark=True)


def modules_mc(out, mc, seed, names):
    weights = np.maximum(mc, 0).astype("float64")
    np.fill_diagonal(weights, 0)
    ids, modularity = bct.community_louvain(weights, gamma=1, seed=seed)
    order = np.argsort(ids, kind="stable")
    # Display every fifth edge; full order and assignments are exported.
    mini = order[::5]
    heatmap(out, "modules_mc", mc[np.ix_(mini,mini)], "Méta-connectivité triée par module", "Lien trié (× 5)", style="mc_m1")
    write_table(out / "data/modules.csv", ["edge_0","roi_i","roi_j","module"], (
        (i,names[EDGE[0][i]],names[EDGE[1][i]],int(ids[i])) for i in order))
    write_json(out / "data/metrics.json", {"Q":float(modularity), "n_modules":len(np.unique(ids))})

