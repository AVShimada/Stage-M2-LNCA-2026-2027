"""Exporter la démographie et les tests cognitifs des sujets avec un second scan."""

import argparse
from datetime import datetime, timezone

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

import commun as io
from cognition import read_rows
from metadonnees import clean_value, column_name, read_demography, session_id, subject_id


def write_workbook(path, name, columns, rows):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = name
    sheet.append(columns)
    for row in rows:
        sheet.append([row[key] for key in columns])
        for cell in sheet[sheet.max_row]:
            if cell.data_type == "f":
                cell.data_type = "s"
    sheet.freeze_panes = "B2"
    sheet.auto_filter.ref = sheet.dimensions
    sheet.sheet_view.showGridLines = False
    for cell in sheet[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="25445E")
        sheet.column_dimensions[get_column_letter(cell.column)].width = max(14, len(cell.value) + 3)
    workbook.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    args = parser.parse_args()

    config = io.ROOT / "config"
    demo_cfg = io.read_json(config / "suivi.json")
    cog_cfg = io.read_json(config / "cognition.json")
    demo_path = (config / demo_cfg["table"]).resolve()
    cog_path = (config / cog_cfg["table"]).resolve()

    demo_columns, demo_rows = read_demography(demo_path, demo_cfg["delimiter"])
    demo_id = column_name(demo_cfg["colonne_sujet"], demo_columns)
    visit2 = demo_cfg["visites"]["ses-2"]
    fc = column_name(visit2["fc_declaree"], demo_columns)
    sc = column_name(visit2["sc_declaree"], demo_columns)
    selected = []
    for row in demo_rows:
        fc_present = clean_value(row[fc], "fc_declaree", None)
        sc_present = clean_value(row[sc], "sc_declaree", None)
        if fc_present is True or sc_present is True:
            selected.append({**row, demo_id: subject_id(row[demo_id])})
    subjects = {row[demo_id] for row in selected}
    for session in demo_cfg["visites"].values():
        for field in ("age", "education"):
            key = column_name(session.get(field), demo_columns)
            if key:
                for row in selected:
                    row[key] = clean_value(row[key], field, None)

    headers, cog_rows = read_rows(cog_path, cog_cfg)
    cog_id = cog_cfg["id_column"]
    visit = cog_cfg["visit_column"]
    cog_columns = [cog_id, visit, cog_cfg["age_column"], "Sex", cog_cfg["education_column"], *cog_cfg["tests"]]
    missing = [key for key in cog_columns if headers.count(key) != 1]
    if missing:
        raise ValueError(f"Colonnes cognitives absentes ou ambiguës : {missing}")
    cognition = []
    for row in cog_rows:
        subject = subject_id(row[cog_id])
        if subject in subjects:
            cognition.append({**row, cog_id: subject})

    cog_subjects = {row[cog_id] for row in cognition}
    counts = {
        "n_sujets_second_scan": len(subjects),
        "n_lignes_demographie": len(selected),
        "n_sujets_avec_cognition": len(cog_subjects),
        "n_sujets_sans_cognition": len(subjects - cog_subjects),
        "visites_cognitives": sorted({session_id(row[visit]) for row in cognition}),
    }
    out = io.stage_dir("cohorte_second_scan") / io.token(args.run)
    title = f"Sujets avec {fc} ou {sc} déclaré, et tests cognitifs disponibles"
    with io.record(out, title, {"fc_seconde_visite": fc, "sc_seconde_visite": sc},
                   {"demographie": demo_path, "cognition": cog_path}, extra=counts) as data:
        write_workbook(data / "second_scan_demographie.xlsx", "Démographie", demo_columns, selected)
        write_workbook(data / "second_scan_cognition.xlsx", "Tests cognitifs", cog_columns, cognition)

    print(f"{len(subjects)} sujets avec second scan, dont {len(cog_subjects)} avec une ligne cognitive")
    if counts["visites_cognitives"]:
        print("Visites cognitives disponibles :", ", ".join(counts["visites_cognitives"]))
    print(out / "data")


if __name__ == "__main__":
    main()
