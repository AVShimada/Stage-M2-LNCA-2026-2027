"""Décrire les sujets avec un second scan à partir des exports du projet."""

import argparse
from collections import Counter
from datetime import datetime, timezone

import numpy as np
from openpyxl import load_workbook

import commun as io
from cohorte import index_table, optional_number
from metadonnees import VISIT_FIELDS, normalize, session_id


BLUE, ORANGE, GREEN = "#4c72b0", "#dd8452", "#55a868"
UNKNOWN = "Non renseigné"
EDUCATION = ("Faible (≤5)", "Moyen (6–8)", "Élevé (>8)")


def completed_export(run, population=None):
    stages = {"filtree": "cohorte_tests_complets", "brute": "cohorte_second_scan"}
    bases = [io.stage_dir(stages[population])] if population else [io.stage_dir(v) for v in stages.values()]
    folders = [base / io.token(run) for base in bases] if run else (p for base in bases for p in base.glob("*"))
    complete = []
    for folder in folders:
        manifest = folder / "manifest.json"
        if manifest.is_file():
            state = io.read_json(manifest)
            if state["status"] == "complete":
                complete.append((folder.parent == io.stage_dir("cohorte_tests_complets"), state["started"], folder))
    if not complete:
        raise ValueError("Aucun export terminé. Exécuter script/extraire_second_scan.py.")
    return max(complete)[-1] / "data/second_scan_demographie.xlsx"


def read_population(path, config):
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        values = iter(workbook["Démographie"].values)
        columns = list(next(values))
        rows = [dict(zip(columns, row)) for row in values if any(value is not None for value in row)]
    finally:
        workbook.close()
    demographic, _, conflicts, n_subjects = normalize(columns, rows, config)
    if not n_subjects:
        raise ValueError("Le classeur démographique ne contient aucun sujet avec un second scan.")
    visits = {"ses-1": {}, "ses-2": {}}
    for subject, session, sex, *fields in demographic:
        visits[session][subject] = dict(zip(VISIT_FIELDS, fields), sex=sex)
    return visits, len(rows), len(conflicts)


def describe(label, values):
    valid = values[np.isfinite(values)]
    if not len(valid):
        return label, 0, len(values), *([None] * 7)
    q1, median, q3 = np.percentile(valid, (25, 50, 75))
    std = valid.std(ddof=1) if len(valid) > 1 else None
    return label, len(valid), len(values) - len(valid), valid.mean(), std, median, q1, q3, valid.min(), valid.max()


def histogram(out, name, series, title, xlabel, bins=30):
    available = [(label, values[np.isfinite(values)], color) for label, values, color in series]
    available = [(label, values, color) for label, values, color in available if len(values)]
    if not available:
        return
    edges = np.histogram_bin_edges(np.concatenate([values for _, values, _ in available]), bins=bins)

    def draw(fig, ax):
        for label, values, color in available:
            ax.hist(values, bins=edges, histtype="stepfilled", alpha=0.6,
                    color=color, label=f"{label} (n = {len(values)})")
        ax.set(title=title, xlabel=xlabel, ylabel="Nombre de sujets")
        ax.legend(frameon=False)

    io.figure(out, name, draw)


def bars(out, name, counts, title, xlabel, color):
    def draw(fig, ax):
        labels = ["Non\nrenseigné" if key == UNKNOWN else key for key in counts]
        chart = ax.bar(labels, list(counts.values()), color=color)
        ax.bar_label(chart, padding=3)
        ax.set(title=title, xlabel=xlabel, ylabel="Nombre de sujets")
        ax.margins(y=0.14)

    io.figure(out, name, draw)


def scan_type(record):
    fc, sc = record["fc_declaree"], record["sc_declaree"]
    if fc is True and sc is True:
        return "FC et SC déclarées"
    if fc is True:
        return "FC déclarée"
    if sc is True:
        return "SC déclarée"
    return UNKNOWN


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohorte-run", help="Identifiant de l’export, dernier export terminé par défaut")
    parser.add_argument("--scores-run", required=True, help="Identifiant du calcul de cognition.py sur le classeur complet")
    parser.add_argument("--run", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    parser.add_argument("--population", choices=["brute", "filtree"], help="Par défaut, priorité à l’export filtré")
    args = parser.parse_args()

    config = io.ROOT / "config"
    demo_path = completed_export(args.cohorte_run, args.population)
    scores_dir = io.stage_dir("cohorte_cognition") / io.token(args.scores_run)
    scores_path = scores_dir / "data/scores.csv"
    score_manifest = io.read_json(scores_dir / "manifest.json")
    cog_cfg = io.read_json(config / "cognition.json")
    cog_source = (config / cog_cfg["table"]).resolve()
    if (score_manifest["status"] != "complete" or score_manifest["configuration"] != cog_cfg
            or score_manifest["inputs"]["table"]["sha256"] != io.digest(cog_source)):
        raise ValueError("Relancer cognition.py sur le classeur source et avec config/cognition.json.")
    if session_id(cog_cfg["session"]) != "ses-1":
        raise ValueError("Le descriptif cognitif utilise les scores de la visite 1.")
    visits, n_rows, n_conflicts = read_population(demo_path, io.read_json(config / "suivi.json"))
    scores = index_table(scores_path, ["cognition"])
    first, second = visits["ses-1"], visits["ses-2"]
    n = len(first)
    age1 = np.array([first[subject]["age"] for subject in first], dtype=float)
    age2 = np.array([second[subject]["age"] for subject in first], dtype=float)
    edu = np.array([first[subject]["education"] for subject in first], dtype=float)
    score = np.array([optional_number(scores.get((subject, "ses-1"), {}), "cognition") for subject in first])
    gap = age2 - age1
    negative = gap < 0
    gap[negative] = np.nan

    threshold = cog_cfg["akt_min_age"]
    age_labels = (f"<{threshold} ans", f"≥{threshold} ans")
    age_group = np.where(np.isnan(age1), UNKNOWN, np.where(age1 < threshold, *age_labels))
    edu_group = np.where(np.isnan(edu), UNKNOWN, np.array(EDUCATION)[np.digitize(edu, (5, 8), right=True)])
    counts = {
        "sexe": Counter(row["sex"] or UNKNOWN for row in first.values()),
        "groupe_age_scan_1": Counter(age_group),
        "ISCED_scan_1": Counter(f"{value:g}" if np.isfinite(value) else UNKNOWN for value in np.sort(edu)),
        "groupe_ISCED_scan_1": Counter(edu_group),
        "modalite_scan_2": Counter(scan_type(row) for row in second.values()),
    }
    numeric = {"age_scan_1": age1, "age_scan_2": age2, "ecart_ages_approx_ans": gap, "score_cognitif_visite_1": score}
    provenance = {
        "n_sujets": n, "n_lignes_demographie": n_rows, "n_champs_en_conflit": n_conflicts,
        "n_ecarts_age_negatifs_exclus": int(negative.sum()), "n_scores_cognitifs": int(np.isfinite(score).sum()),
    }
    out = io.stage_dir("cohorte_second_scan_descriptif") / io.token(args.run)
    settings = {"seuil_age_m1": threshold, "seuils_isced_m1": [5, 8], "scores_run": args.scores_run,
                "delai": "Différence des âges arrondis entre les scans, en années"}
    with io.record(out, "Descriptifs des sujets avec un second scan", settings,
                   {"demographie": demo_path, "scores_visite_1": scores_path,
                    "manifest_cognition": scores_dir / "manifest.json"}, extra=provenance) as data:
        io.write_table(data / "statistiques.csv", ["variable", "n", "n_manquant", "moyenne", "ecart_type",
                       "mediane", "q1", "q3", "minimum", "maximum"],
                       (describe(label, values) for label, values in numeric.items()))
        io.write_table(data / "effectifs.csv", ["variable", "categorie", "n", "pct_sur_total"],
                       ((key, label, count, 100 * count / n) for key, groups in counts.items() for label, count in groups.items()))
        histogram(out, "age_scans", [("Scan 1", age1, BLUE), ("Scan 2", age2, ORANGE)],
                  "Distribution des âges aux deux scans", "Âge (années)")
        histogram(out, "age_groupes_scan_1", [(label, age1[age_group == label], color)
                  for label, color in zip(age_labels, (BLUE, ORANGE))],
                  "Âge au scan 1 par groupe", "Âge (années)")
        bars(out, "sexe", counts["sexe"], "Répartition par sexe", "Sexe déclaré", BLUE)
        bars(out, "education", counts["ISCED_scan_1"], "Niveau d’éducation au scan 1", "ISCED", GREEN)
        histogram(out, "delai_approx", [("Âges disponibles", gap, BLUE)],
                  "Délai approximatif entre les deux scans", "Âge scan 2 − âge scan 1 (années)")
        histogram(out, "score_cognition", [("Visite 1", score, BLUE)],
                  "Distribution du score global de cognition", "Score global")

        def cognition_age(fig, ax):
            for group, color in zip(age_labels, (BLUE, ORANGE)):
                keep = (age_group == group) & np.isfinite(score)
                if keep.any():
                    ax.scatter(age1[keep], score[keep], color=color, alpha=0.65, label=f"{group} (n = {keep.sum()})")
            ax.set(title="Score global de cognition selon l’âge", xlabel="Âge au scan 1 (années)", ylabel="Score global · visite 1")
            ax.legend(frameon=False)

        if (np.isfinite(age1) & np.isfinite(score)).any():
            io.figure(out, "cognition_age", cognition_age)

        def cognition_education(fig, ax):
            groups = [(label, score[(edu_group == label) & np.isfinite(score)]) for label in EDUCATION]
            groups = [(label, values) for label, values in groups if len(values)]
            ax.boxplot([values for _, values in groups], patch_artist=True, boxprops={"facecolor": "#a6bde0"})
            ax.set_xticks(range(1, len(groups) + 1), [f"{label}\nn = {len(values)}" for label, values in groups])
            ax.set(title="Score global de cognition selon le niveau d’éducation", xlabel="Groupe ISCED au scan 1",
                   ylabel="Score global · visite 1")

        if (np.isfinite(edu) & np.isfinite(score)).any():
            io.figure(out, "cognition_education", cognition_education)

    print(f"{n} sujets avec un second scan, {provenance['n_scores_cognitifs']} scores cognitifs disponibles")
    print(out)


if __name__ == "__main__":
    main()
