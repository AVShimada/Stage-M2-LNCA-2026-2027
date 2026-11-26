"""Score global M1 et ACP descriptive sur la table cognitive d'une visite."""
import argparse
import csv
from pathlib import Path

import numpy as np
import commun as io
from metadonnees import session_id, subject_id


def read_rows(path, config):
    if path.suffix.lower() != ".xlsx":
        io.required(config, ["delimiter"])
        with path.open(encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f, delimiter=config["delimiter"])
            return reader.fieldnames or [], list(reader)

    from openpyxl import load_workbook

    io.required(config, ["sheet", "visit_column"])
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        values = iter(workbook[config["sheet"]].values)
        headers = [str(value).strip() if value is not None else "" for value in next(values)]
        return headers, [dict(zip(headers, row)) for row in values if any(value is not None for value in row)]
    finally:
        workbook.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",required=True,type=Path)
    p.add_argument("--run",required=True)
    args=p.parse_args()
    c=io.read_json(args.config)
    io.required(c,["table","id_column","age_column","education_column","min_tests","session"])
    src=(args.config.parent / c["table"]).resolve()
    keys=list(c["tests"])
    headers, rows = read_rows(src, c)
    visit = c.get("visit_column")
    columns = [c[k] for k in ("id_column", "age_column", "education_column")] + keys
    if visit:
        columns.append(visit)
    missing_columns = [name for name in columns if name not in headers]
    if missing_columns:
        raise ValueError(f"Colonnes cognitives absentes : {missing_columns}")
    if any(headers.count(name) != 1 for name in columns):
        raise ValueError("Colonnes cognitives en double.")
    session = session_id(c["session"])
    if visit:
        rows = [row for row in rows if session_id(row[visit]) == session]
    ids = np.array([subject_id(row[c["id_column"]]) for row in rows])
    if not rows or len(ids) != len(set(ids)):
        raise ValueError("Fournir une ligne par sujet pour la session indiquée.")
    if not 1 <= c["min_tests"] <= len(keys) or not set(c["tests"].values()) <= {-1,1}:
        raise ValueError("Vérifier le minimum de tests et leur orientation")
    missing = {"", *(str(v).strip().casefold() for v in c["missing_codes"])}
    def value(r,k):
        text = str(r[k]).strip() if r[k] is not None else ""
        return float(text.replace(",", ".")) if text.casefold() not in missing else np.nan
    age=np.array([value(r,c["age_column"]) for r in rows])
    edu=np.array([value(r,c["education_column"]) for r in rows])
    x=np.array([[value(r,k) for k in keys] for r in rows])
    missing_age = np.isnan(age)
    if c.get("drop_missing_age"):
        ids, age, edu, x = (a[~missing_age] for a in (ids, age, edu, x))
    if not len(ids) or not np.isfinite(age).all() or np.isinf(edu).any() or np.isinf(x).any():
        raise ValueError("Vérifier âge, éducation et tests non finis")
    if "AKT_TRW_SelectiveAttention" in keys:
        x[age<c["akt_min_age"],keys.index("AKT_TRW_SelectiveAttention")]=np.nan
    counts=np.isfinite(x).sum(axis=1)
    if np.any(np.sum(np.isfinite(x),axis=0)<2) or not np.all(np.nanstd(x,axis=0)>0):
        raise ValueError("Test absent ou constant sur la table de référence")
    mean,std=np.nanmean(x,axis=0),np.nanstd(x,axis=0)
    signed=(x-mean)/std*np.array(list(c["tests"].values()))
    score=np.divide(np.nansum(signed,axis=1),counts,out=np.full(len(x),np.nan),where=counts>=c["min_tests"])
    # ACP descriptive avec imputation par moyenne puis standardisation, comme au M1.
    imputed=np.where(np.isnan(x),mean,x)
    zz=(imputed-imputed.mean(axis=0))/imputed.std(axis=0)
    u,s,vt=np.linalg.svd(zz,full_matrices=False)
    variance=s**2/np.sum(s**2)
    folder=io.stage_dir("cohorte_cognition")/io.token(args.run)
    provenance = {"n_lignes_visite": len(rows), "n_age_manquant_exclus": len(rows) - len(ids),
                  "n_scores": int(np.isfinite(score).sum())}
    with io.record(folder,"Score cognitif global et ACP descriptive",c,{"table":src}, extra=provenance) as data:
        io.write_table(data/"scores.csv",["subject","session","age","education","cognition","tests_used"],
                       zip(ids,[session]*len(ids),age,edu,score,counts))
        io.write_table(data/"standardisation.csv",["test","sign","mean","std","n"],zip(keys,c["tests"].values(),mean,std,np.isfinite(x).sum(axis=0)))
        pcs=[f"PC{i+1}" for i in range(len(s))]
        io.write_table(data/"acp_scores.csv",["subject",*pcs],([i,*v] for i,v in zip(ids,u*s)))
        io.write_table(data/"acp_axes.csv",["test",*pcs],([k,*v] for k,v in zip(keys,vt.T)))
        io.write_table(data/"acp_variance.csv",["component","fraction"],zip(pcs,variance))
        def draw(fig,ax):
            ax.plot(np.arange(1,len(s)+1),np.cumsum(variance),"o-",color="#4c72b0")
            ax.set(xlabel="Nombre de composantes",ylabel="Fraction de variance cumulée",title="ACP cognitive descriptive")
        io.figure(folder,"acp",draw)
    print(f"{len(ids)} sujets retenus en {session}, {provenance['n_scores']} scores calculés")
    print(folder)


if __name__=="__main__":
    main()

