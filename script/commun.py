"""Entrées et traçabilité communes au pipeline 1000BRAINS."""
import csv
import hashlib
import json
import os
import platform
import re
import shlex
import sys
import traceback
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
EDGE = np.triu_indices(100, 1)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def write_table(path, header, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def token(value):
    if not isinstance(value, str) or not re.fullmatch(r"[\w.-]+", value) or value in (".", ".."):
        raise ValueError(f"Identifiant vide ou invalide : {value!r}")
    return value


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def required(cfg, keys):
    missing = [key for key in keys if cfg.get(key) in (None, "")]
    if missing:
        raise ValueError("À compléter : " + ", ".join(missing))


def source(cfg, key):
    required(cfg, [key])
    p = Path(cfg[key]).expanduser()
    if not p.is_absolute():
        required(cfg, ["dataset_root"])
        base = Path(cfg["dataset_root"]).expanduser()
        if not base.is_absolute():
            raise ValueError("dataset_root doit être absolu")
        p = base / p
    if not p.is_file():
        raise FileNotFoundError(f"{key} : {p}")
    return p.resolve()


def array_file(path):
    if path.suffix.lower() == ".npy":
        return np.load(path, allow_pickle=False)
    delim = "," if path.suffix.lower() == ".csv" else None
    return np.loadtxt(path, delimiter=delim, encoding="utf-8-sig")


def labels(path):
    names = [s.strip() for s in path.read_text(encoding="utf-8-sig").splitlines()]
    if len(names) != 100 or len(set(names)) != 100 or "" in names:
        raise ValueError(f"Il faut 100 labels distincts, un par ligne : {path}")
    return names


def aligned_sc(path, names_path, reference):
    names = labels(names_path)
    if set(names) != set(reference):
        raise ValueError("Les régions SC et BOLD diffèrent")
    index = [names.index(name) for name in reference]
    sc = np.asarray(array_file(path), dtype=float)
    if sc.shape != (100, 100) or not np.isfinite(sc).all() or not np.allclose(sc, sc.T, atol=1e-8):
        raise ValueError("SC attendue : 100 × 100, finie et symétrique")
    sc = sc[np.ix_(index, index)]
    np.fill_diagonal(sc, 0)
    return sc


def load_inputs(c):
    required(c, ["subject", "session", "version", "tr_seconds", "bold", "sc", "labels_bold",
                 "labels_sc", "sc_weights", "bold_orientation", "censor_mode"])
    if c["atlas"] != "Schaefer 100" or c["version"] not in ("V1", "V4"):
        raise ValueError("Atlas ou version inattendus")
    tr = float(c["tr_seconds"])
    if not np.isfinite(tr) or tr <= 0:
        raise ValueError("TR non valide")
    files = {k: source(c, k) for k in ("bold", "sc", "labels_bold", "labels_sc")}
    for mark in re.findall(r"(?:^|[_\-.])(V[1-9]\d*)(?=[_\-.]|$)", files["bold"].name, re.I):
        if mark.upper() != c["version"]:
            raise ValueError(f"Version BOLD contradictoire : {files['bold'].name}")
    names = labels(files["labels_bold"])
    sc = aligned_sc(files["sc"], files["labels_sc"], names)
    bold = np.asarray(array_file(files["bold"]), dtype=float)
    if c["bold_orientation"] == "regions_time":
        bold = bold.T
    elif c["bold_orientation"] != "time_regions":
        raise ValueError("bold_orientation : time_regions ou regions_time")
    if bold.ndim != 2 or bold.shape[1] != 100:
        raise ValueError(f"Format BOLD attendu : temps × 100 après orientation, reçu {bold.shape}")
    mode = c["censor_mode"]
    n = len(bold)
    grid = np.arange(n)
    if mode in ("keep_mask", "retained_indices"):
        files["censor_file"] = source(c, "censor_file")
        v = np.atleast_1d(array_file(files["censor_file"]))
        if v.ndim != 1 or not np.isfinite(v).all() or not np.all(v == np.floor(v)):
            raise ValueError("Censure attendue sous forme d'un vecteur entier")
        if mode == "keep_mask":
            if len(v) != n or not np.isin(v, [0, 1]).all():
                raise ValueError("Masque requis sur toute la série : 1 = garder, 0 = exclure")
            bold, grid = bold[v == 1], grid[v == 1]
        else:
            required(c, ["original_volume_count"])
            n = c["original_volume_count"]
            grid = v.astype(int)
            if type(n) is not int or len(grid) != len(bold) or np.any(np.diff(grid) <= 0) or np.any((grid < 0) | (grid >= n)):
                raise ValueError("Indices requis : croissants, base zéro, un par volume BOLD conservé")
    elif mode != "none":
        raise ValueError("censor_mode : none, keep_mask ou retained_indices")
    if len(bold) < 3 or not np.isfinite(bold).all() or np.any(np.ptp(bold, axis=0) == 0):
        raise ValueError("Volumes retenus insuffisants, non finis ou région constante")
    return bold, grid, n, sc, names, files


def run_dir(c, stage, run):
    tag = "_".join(token(x) for x in (run, c["subject"], c["session"], c["version"]))
    return OUTPUT / stage / tag


def code_hashes():
    return {p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def versions():
    result = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "numba", "bctpy", "networkx", "matplotlib"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            pass
    return result


@contextmanager
def record(out, title, cfg, inputs, extra=None):
    out.parent.mkdir(parents=True, exist_ok=True)
    description = out.parent / "README.md"
    try:
        with description.open("x", encoding="utf-8") as f:
            f.write(f"# {out.parent.name}\n\n{title}\n\nVoir WORKFLOW_1000BRAINS.md pour la méthode et les prérequis. "
                    f"Chaque exécution conserve ses entrées et ses paramètres dans manifest.json, ses valeurs dans data/ et ses figures dans figures/.\n"
                    f"\nExemple du premier lancement : `{shlex.join([sys.executable,*sys.argv])}`\n")
    except FileExistsError:
        pass
    out.mkdir(exist_ok=False)
    for folder in ("figures", "data", "logs"):
        (out / folder).mkdir()
    m = dict(status="running", started=stamp(), configuration=cfg, command=sys.argv, code_sha256=code_hashes(),
             environment=versions(), inputs={k: {"path": str(v), "sha256": digest(v)} for k, v in inputs.items()},
             job_id=next((os.getenv(k) for k in ("SLURM_JOB_ID", "PBS_JOBID", "JOB_ID", "LSB_JOBID") if os.getenv(k)), None))
    if extra:
        m.update(extra)
    write_json(out / "manifest.json", m)
    with (out / "logs/run.log").open("w", encoding="utf-8") as log:
        try:
            with redirect_stdout(log), redirect_stderr(log):
                print(title, flush=True)
                yield out / "data"
            m["status"] = "complete"
        except BaseException:
            m["status"] = "failed"
            traceback.print_exc(file=log)
            raise
        finally:
            m["finished"] = stamp()
            write_json(out / "manifest.json", m)
            (out / "README.md").write_text(f"# {title}\n\nStatut : {m['status']}\n\nCommande : `{ ' '.join(sys.argv) }`\n\nDétails dans manifest.json et logs/run.log.\n", encoding="utf-8")


def completed(c, stage, run):
    out = run_dir(c, stage, run)
    m = read_json(out / "manifest.json")
    saved, current = m["configuration"].copy(), c.copy()
    if "coordinates_csv" not in m["inputs"]:
        saved.pop("coordinates_csv", None)
        current.pop("coordinates_csv", None)
    if m["status"] != "complete" or saved != current or m["code_sha256"] != code_hashes():
        raise ValueError(f"Exécution incompatible ou incomplète : {out}")
    return out, m


MODEL_COLORS = {"mfm": "#0072BD", "tmfm": "#D95319"}
SPEED_STYLES = {"emp": ((0.1, 0.1, 0.5), "-"), "phase": ("g", "--"), "shuffle": ("r", ":")}


def figure(out, name, draw, *, dark=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.style.context("dark_background" if dark else "default"):
        fig, ax = plt.subplots(figsize=(7.6, 5.2), layout="constrained")
        try:
            draw(fig, ax)
            for ext in ("png", "pdf"):
                fig.savefig(out / "figures" / f"{name}.{ext}", dpi=150)
        finally:
            plt.close(fig)


def heatmap(out, name, matrix, title, axes="Région", signed=True, *, style=None, limits=None):
    # Le bloc historique du pilote appelle simplement heatmap(..., "fcd", ...).
    if style is None:
        style = "fcd_article" if name == "fcd" else "matrice_m1"
    palettes = {
        "matrice_m1": ("turbo", (-1, 1) if signed else (None, None), False),
        "fcd_m1": ("turbo", (0, 0.5), True),
        "fcd_article": ("jet", (0, 1), False),
        "mc_m1": ("turbo", (0, 0.1), True),
        "transitions_fc": ("Blues", (0, None), False),
        "transitions_sc": ("Greens", (0, None), False),
    }
    cmap, bounds, dark = palettes[style]
    low, high = bounds if limits is None else limits
    origin = "upper" if style in ("fcd_m1", "fcd_article") else "lower"
    def draw(fig, ax):
        im = ax.imshow(matrix, aspect="auto", interpolation="nearest", origin=origin,
                       cmap=cmap, vmin=low, vmax=high)
        ax.set(xlabel=axes, ylabel=axes, title=title)
        below = low is not None and np.nanmin(matrix) < low
        above = high is not None and np.nanmax(matrix) > high
        extend = {(False,False):"neither", (True,False):"min",
                  (False,True):"max", (True,True):"both"}[below, above]
        fig.colorbar(im, ax=ax, extend=extend)
    figure(out, name, draw, dark=dark)
