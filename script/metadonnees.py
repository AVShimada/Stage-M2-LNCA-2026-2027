"""Importer la démographie 1000BRAINS et inventorier le fichier MATLAB."""
import argparse
import csv
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path
import re

import numpy as np

import commun as io

# Codes de valeur indisponible dans les exports Excel et CSV.
MISSING = {"", "na", "n/a", "nan", "null", "#nv", "#n/a", "#n/d"}
VISIT_FIELDS = ("age", "education", "fc_declaree", "sc_declaree", "date_scan")


def subject_id(value):
    text = str(value).strip() if value is not None else ""
    match = re.fullmatch(r"(?:sub-)?([0-9]+)(?:\.0+)?", text, re.IGNORECASE)
    if match:
        return f"sub-{int(match[1]):04d}"
    if text.casefold() in MISSING:
        raise ValueError("Identifiant sujet manquant.")
    return io.token(text)


def session_id(value):
    text = str(value).strip() if value is not None else ""
    match = re.fullmatch(r"(?:ses-)?([0-9]+)(?:\.0+)?", text, re.IGNORECASE)
    if match:
        return f"ses-{int(match[1])}"
    if text.casefold() in MISSING:
        raise ValueError("Identifiant de visite manquant.")
    return io.token(text)


def read_demography(path, delimiter):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"En-tête démographique absent : {path}")
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Colonnes démographiques en double.")
        return reader.fieldnames, list(reader)


def column_name(name, headers):
    if not name or name in headers:
        return name
    compact = lambda s: re.sub(r"[^a-z0-9]", "", s.casefold())
    matches = [header for header in headers if compact(header) == compact(name)]
    if len(matches) != 1:
        raise ValueError(f"Colonne démographique absente ou ambiguë : {name}. En-têtes : {headers}")
    return matches[0]


def clean_value(value, field, date_format):
    text = str(value).strip() if value is not None else ""
    if text.casefold() in MISSING:
        return None
    if field in ("age", "education", "intervalle"):
        number = float(text.replace(",", "."))
        if not np.isfinite(number):
            raise ValueError(f"Valeur non finie pour {field}.")
        return number
    if field in ("fc_declaree", "sc_declaree"):
        flags = {"yes": True, "no": False, "oui": True, "non": False,
                 "1": True, "0": False, "true": True, "false": False}
        try:
            return flags[text.casefold()]
        except KeyError:
            raise ValueError(f"Indicateur de présence inconnu pour {field}.") from None
    if field == "date_scan":
        parsed = datetime.strptime(text, date_format).date() if date_format else date.fromisoformat(text)
        return parsed.isoformat()
    return text


def mat_inventory(path):
    from scipy.io import loadmat, whosmat

    try:
        return {"variables": [
            {"nom": name, "dimensions": list(shape), "classe": kind}
            for name, shape, kind in whosmat(path)
        ]}
    except TypeError:
        # Certains objets MATLAB n'ont pas de dimensions dans leur en-tête.
        values = loadmat(path)
        return {"methode": "loadmat", "variables": [
            {"nom": name, "dimensions": list(getattr(value, "shape", ())),
             "classe": type(value).__name__}
            for name, value in values.items() if name and not name.startswith("__")
        ]}
    except NotImplementedError:
        try:
            import h5py
        except ImportError:
            return {"format": "MATLAB v7.3", "inspection": "Installer h5py pour inventorier les variables."}
        with h5py.File(path) as mat:
            return {"format": "MATLAB v7.3", "variables": [
                {"nom": name, "dimensions": list(value.shape) if isinstance(value, h5py.Dataset) else [],
                 "classe": "dataset" if isinstance(value, h5py.Dataset) else "group"}
                for name, value in mat.items() if not name.startswith("#")
            ]}


def normalize(headers, rows, cfg):
    key = column_name(cfg["colonne_sujet"], headers)
    sex = column_name(cfg["colonne_sexe"], headers)
    interval = column_name(cfg.get("colonne_intervalle"), headers)
    visits = ("ses-1", "ses-2")
    mapping = {
        session: {field: column_name(cfg["visites"][session].get(field), headers) for field in VISIT_FIELDS}
        for session in visits
    }
    factors = {"jours": 1 / 365.25, "mois": 1 / 12, "ans": 1}
    factor = factors[cfg["unite_intervalle"]] if interval else None
    grouped = defaultdict(list)
    for row in rows:
        grouped[subject_id(row[key])].append(row)
    demographic, followup, conflicts = [], [], []
    for subject, group in sorted(grouped.items()):
        def value(column, field, session=""):
            if not column:
                return None
            found = {clean_value(row[column], field, cfg.get("format_date")) for row in group}
            found.discard(None)
            if len(found) > 1:
                conflicts.append((subject, session, field))
                return None
            return next(iter(found), None)

        gender = value(sex, "sex")
        by_visit = {}
        for session in visits:
            values = {field: value(column, field, session) for field, column in mapping[session].items()}
            by_visit[session] = values
            demographic.append((subject, session, gender, *values.values()))
        first, second = (by_visit[session] for session in visits)
        start, end = first["date_scan"], second["date_scan"]
        years, origin = None, ""
        if start and end:
            days = (date.fromisoformat(end) - date.fromisoformat(start)).days
            if days > 0:
                years, origin = days / 365.25, "dates_scan"
            else:
                conflicts.append((subject, "", "ordre_dates_scan"))
        elif interval:
            length = value(interval, "intervalle")
            if length is not None and length > 0:
                years, origin = length * factor, interval
            elif length is not None:
                conflicts.append((subject, "", "intervalle_non_positif"))
        followup.append((subject, first["fc_declaree"], first["sc_declaree"], second["fc_declaree"],
                         second["sc_declaree"], start, end, years, origin))
    return demographic, followup, conflicts, len(grouped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/suivi.json"))
    parser.add_argument("--run", required=True)
    args = parser.parse_args()
    cfg = io.read_json(args.config)
    table, mat = ((args.config.parent / cfg[key]).resolve() for key in ("table", "scores_mat"))
    headers, rows = read_demography(table, cfg["delimiter"])
    demographic, followup, conflicts, subjects = normalize(headers, rows, cfg)
    if not subjects:
        raise ValueError("La table démographique ne contient aucun sujet.")
    # L'inventaire MATLAB est indépendant des informations démographiques.
    from scipy.io.matlab import MatReadError

    try:
        inventory = mat_inventory(mat)
    except (MatReadError, TypeError, ValueError, OSError, NotImplementedError) as error:
        inventory = {"variables": [], "inspection": f"{type(error).__name__}: {error}"}
    summary = {"n_lignes_source": len(rows), "n_sujets": subjects, "n_conflits": len(conflicts),
               "n_intervalles": sum(row[-2] is not None for row in followup), "colonnes_source": headers}
    out = io.stage_dir("cohorte_metadonnees") / io.token(args.run)
    with io.record(out, "Métadonnées démographiques et visites 1000BRAINS", cfg,
                   {"demographie": table, "scores_mat": mat}) as data:
        io.write_table(data / "demographie.csv", ["subject", "session", "sex", *VISIT_FIELDS], demographic)
        io.write_table(data / "suivi.csv", ["subject", "fc_ses1", "sc_ses1", "fc_ses2", "sc_ses2",
                       "date_ses1", "date_ses2", "intervalle_annees", "source_intervalle"], followup)
        io.write_table(data / "conflits.csv", ["subject", "session", "champ"], conflicts)
        io.write_json(data / "resume.json", summary)
        io.write_json(data / "cog_scores_mat.json", inventory)
    print(f"{subjects} sujets, {len(conflicts)} champs à vérifier, {summary['n_intervalles']} intervalles renseignés")
    if "inspection" in inventory:
        print(f"Inventaire MATLAB incomplet : {inventory['inspection']}")
    print(out)


if __name__ == "__main__":
    main()
