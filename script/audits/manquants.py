"""Compter les données manquantes et les scores à zéro dans les métadonnées."""

import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import commun as io
from cognition import read_rows
from metadonnees import MISSING, clean_value, column_name, read_demography, session_id, subject_id


def code(value):
    return str(value).strip().casefold() if value is not None else ""


def number(value, missing):
    text = code(value)
    if text in missing:
        return np.nan, False
    try:
        parsed = float(text.replace(",", "."))
    except ValueError:
        return np.nan, True
    return (parsed, False) if np.isfinite(parsed) else (np.nan, True)


def cognitive_audit(rows, cfg, missing):
    tests = list(cfg["tests"])
    columns = [cfg["age_column"], *tests]
    parsed = np.array([[number(row[key], missing) for key in columns] for row in rows], dtype=float)
    parsed = parsed.reshape(len(rows), len(columns), 2)
    age, values = parsed[:, 0, 0], parsed[:, 1:, 0]
    invalid = parsed[:, 1:, 1].astype(bool)
    null = np.isnan(values) & ~invalid
    known_age = np.isfinite(age)
    expected = np.broadcast_to(known_age[:, None], values.shape).copy()
    if "AKT_TRW_SelectiveAttention" in tests:
        expected[:, tests.index("AKT_TRW_SelectiveAttention")] &= age >= cfg["akt_min_age"]
    counts = np.sum(expected & np.isfinite(values), axis=1)[known_age]

    summary = {
        "n_lignes_cognitives": len(rows),
        "n_sujets_cognitifs": len({subject_id(row[cfg["id_column"]]) for row in rows if code(row[cfg["id_column"]]) not in missing}),
        "n_ages_indisponibles": int((~known_age).sum()),
        "n_nulls_par_champ": {key: sum(code(row[key]) in missing for row in rows)
                              for key in (cfg["id_column"], cfg["age_column"], "Sex", cfg["education_column"]) if not rows or key in rows[0]},
        "n_cellules_tests_manquantes": int(null.sum()),
        "n_cellules_tests_non_numeriques_ou_non_finies": int(invalid.sum()),
        "n_cellules_tests_a_zero": int((values == 0).sum()),
        "n_lignes_avec_au_moins_un_test_manquant": int(null.any(axis=1).sum()),
        "n_lignes_tous_tests_manquants": int(null.all(axis=1).sum()),
        "n_lignes_age_connu_avec_test_attendu_indisponible": int((expected & ~np.isfinite(values)).any(axis=1).sum()),
        "n_lignes_age_connu_sous_min_tests": int((counts < cfg["min_tests"]).sum()),
        "distribution_n_tests_utilisables_age_connu": dict(sorted(Counter(counts.tolist()).items())),
        "codes_manquants_tests": dict(Counter(code(row[key]) or "(vide)" for row in rows for key in tests
                                               if code(row[key]) in missing)),
    }
    details = [(name, len(rows), int(null[:, j].sum()), int(invalid[:, j].sum()), int((values[:, j] == 0).sum()),
                int(expected[:, j].sum()), int((expected[:, j] & null[:, j]).sum()),
                int((expected[:, j] & invalid[:, j]).sum())) for j, name in enumerate(tests)]
    return summary, details


def demographic_groups(headers, rows, cfg):
    key = column_name(cfg["colonne_sujet"], headers)
    visit2 = cfg["visites"]["ses-2"]
    fc = column_name(visit2["fc_declaree"], headers)
    sc = column_name(visit2["sc_declaree"], headers)
    indexed = [(subject_id(row[key]) if code(row[key]) not in MISSING else None, row) for row in rows]
    subjects = {
        subject for subject, row in indexed if subject is not None and (
            clean_value(row[fc], "fc_declaree", None) is True or clean_value(row[sc], "sc_declaree", None) is True)
    }
    return subjects, {"ensemble": rows, "second_scan": [row for subject, row in indexed if subject in subjects]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    parser.add_argument("--cognition", type=Path, help="Chemin du classeur si différent de config/cognition.json")
    parser.add_argument("--cognition-seule", action="store_true")
    args = parser.parse_args()

    config = io.ROOT / "config"
    cfg = io.read_json(config / "cognition.json")
    src = (args.cognition or config / cfg["table"]).resolve()
    headers, records = read_rows(src, cfg)
    visit, key = cfg["visit_column"], cfg["id_column"]
    required = [key, visit, cfg["age_column"], cfg["education_column"], *cfg["tests"]]
    invalid_headers = [name for name in required if headers.count(name) != 1]
    if invalid_headers:
        raise ValueError(f"Colonnes cognitives absentes ou ambiguës : {invalid_headers}")
    missing = MISSING | {code(value) for value in cfg["missing_codes"]}
    session = session_id(cfg["session"])
    rows = [row for row in records if code(row[visit]) not in missing and session_id(row[visit]) == session]
    groups = {"ensemble": rows}
    inputs = {"cognition": src}
    report = {"session_cognitive": session, "n_lignes_code_visite_manquant": sum(code(row[visit]) in missing for row in records)}
    demo_details = []
    if not args.cognition_seule:
        demo_cfg = io.read_json(config / "suivi.json")
        demo_src = (config / demo_cfg["table"]).resolve()
        demo_headers, demo_rows = read_demography(demo_src, demo_cfg["delimiter"])
        subjects, demo_groups = demographic_groups(demo_headers, demo_rows, demo_cfg)
        groups["second_scan"] = [row for row in rows if code(row[key]) not in missing and subject_id(row[key]) in subjects]
        demo_details = [(group, name, len(values), sum(code(row[name]) in MISSING for row in values))
                        for group, values in demo_groups.items() for name in demo_headers]
        report["n_sujets_second_scan_declares"] = len(subjects)
        report["n_sujets_second_scan_sans_ligne_cognitive"] = len(subjects - {subject_id(row[key]) for row in groups["second_scan"]})
        inputs["demographie"] = demo_src

    details = []
    for scope, sample in groups.items():
        report[scope], tests = cognitive_audit(sample, cfg, missing)
        details.extend((scope, *test) for test in tests)
    settings = {"session": session, "akt_min_age": cfg["akt_min_age"], "min_tests": cfg["min_tests"]}
    out = io.stage_dir("audit_manquants") / io.token(args.run)
    with io.record(out, "Valeurs manquantes dans les métadonnées", settings, inputs) as data:
        io.write_json(data / "resume.json", report)
        io.write_table(data / "cognition_par_test.csv", ["population", "test", "n_lignes", "n_manquant",
                       "n_invalide", "n_zero", "n_attendu_age_connu", "n_manquant_attendu", "n_invalide_attendu"], details)
        if demo_details:
            io.write_table(data / "demographie_par_champ.csv", ["population", "champ", "n_lignes", "n_manquant"], demo_details)
    for scope in groups:
        stat = report[scope]
        print(f"{scope} : {stat['n_sujets_cognitifs']} sujets avec ligne cognitive, "
              f"{stat['n_lignes_avec_au_moins_un_test_manquant']} avec au moins un test manquant")
    print(out / "data")


if __name__ == "__main__":
    main()
