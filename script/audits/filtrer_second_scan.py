"""Créer les classeurs des sujets suivis dont les tests attendus sont complets."""

import argparse
from collections import Counter
from datetime import datetime, timezone
from math import isfinite

from openpyxl import load_workbook

import commun as io
from extraire_second_scan import write_workbook
from manquants import code, number
from metadonnees import MISSING, column_name, session_id, subject_id


AKT = "AKT_TRW_SelectiveAttention"


def source_run(base, run):
    paths = [base / io.token(run)] if run else base.glob("*")
    completed = []
    for path in paths:
        manifest = path / "manifest.json"
        if manifest.is_file():
            state = io.read_json(manifest)
            if state["status"] == "complete" and not state["configuration"].get("filtre_tests_attendus"):
                completed.append((state["started"], path))
    if not completed:
        raise ValueError(f"Aucun export initial terminé dans {base}.")
    return max(completed)[1]


def read_sheet(path, name):
    book = load_workbook(path, read_only=True, data_only=True)
    try:
        values = iter(book[name].values)
        columns = list(next(values))
        rows = [dict(zip(columns, row)) for row in values if any(value is not None for value in row)]
        return columns, rows
    finally:
        book.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", help="Identifiant de l’export initial, le plus récent par défaut")
    parser.add_argument("--run", help="Identifiant du lancement de sortie")
    args = parser.parse_args()

    base = io.stage_dir("cohorte_second_scan")
    source = source_run(base, args.source_run)
    demo_path = source / "data/second_scan_demographie.xlsx"
    cog_path = source / "data/second_scan_cognition.xlsx"
    demo_cols, demo_rows = read_sheet(demo_path, "Démographie")
    cog_cols, cog_rows = read_sheet(cog_path, "Tests cognitifs")

    config = io.ROOT / "config"
    cfg = io.read_json(config / "cognition.json")
    demo_cfg = io.read_json(config / "suivi.json")
    demo_id = column_name(demo_cfg["colonne_sujet"], demo_cols)
    cog_id, visit, age_col = (cfg[key] for key in ("id_column", "visit_column", "age_column"))
    tests = list(cfg["tests"])
    unavailable = [key for key in (cog_id, visit, age_col, *tests) if cog_cols.count(key) != 1]
    if unavailable:
        raise ValueError(f"Colonnes cognitives absentes ou ambiguës : {unavailable}")
    session = session_id(cfg["session"])
    missing = MISSING | {code(value) for value in cfg["missing_codes"]}
    baseline = {}
    for row in cog_rows:
        if session_id(row[visit]) == session:
            subject = subject_id(row[cog_id])
            if subject in baseline:
                raise ValueError("Plusieurs lignes cognitives pour un même sujet à la visite de référence.")
            baseline[subject] = row

    subjects = {subject_id(row[demo_id]) for row in demo_rows}
    keep = set()
    reasons, missing_tests = Counter(), Counter()
    for subject in subjects:
        row = baseline.get(subject)
        if row is None:
            reasons["sans_ligne_cognitive"] += 1
            continue
        age, _ = number(row[age_col], missing)
        if not isfinite(age):
            reasons["age_cognitif_indisponible"] += 1
            continue
        expected = [test for test in tests if test != AKT or age >= cfg["akt_min_age"]]
        absent = [test for test in expected if not isfinite(number(row[test], missing)[0])]
        if absent:
            reasons["test_attendu_indisponible"] += 1
            missing_tests.update(absent)
        else:
            keep.add(subject)

    kept_demo = [row for row in demo_rows if subject_id(row[demo_id]) in keep]
    kept_cog = [row for row in cog_rows if subject_id(row[cog_id]) in keep]
    run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "_tests_complets"
    out = io.stage_dir("cohorte_tests_complets") / io.token(args.run or run)
    settings = {"source_run": source.name, "filtre_tests_attendus": True,
                "session_cognitive": session, "akt_min_age": cfg["akt_min_age"]}
    summary = {"n_sujets_source": len(subjects), "n_sujets_retenus": len(keep),
               "n_exclus": len(subjects - keep), "exclusions_par_motif": dict(reasons),
               "test_attendu_indisponible_par_test": dict(sorted(missing_tests.items()))}
    with io.record(out, "Sous-groupe second scan avec tous les tests attendus", settings,
                   {"demographie": demo_path, "cognition": cog_path}, extra=summary) as data:
        write_workbook(data / demo_path.name, "Démographie", demo_cols, kept_demo)
        write_workbook(data / cog_path.name, "Tests cognitifs", cog_cols, kept_cog)
        io.write_json(data / "resume_filtrage.json", summary)

    print(f"{len(keep)}/{len(subjects)} sujets retenus, {len(subjects - keep)} exclus")
    print(out / "data")


if __name__ == "__main__":
    main()
