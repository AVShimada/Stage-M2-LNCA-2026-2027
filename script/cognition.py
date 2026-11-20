"""Score global M1 et ACP descriptive sur la table cognitive d'une visite."""
import argparse
import csv
from pathlib import Path

import numpy as np
import commun as io


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",required=True,type=Path)
    p.add_argument("--run",required=True)
    args=p.parse_args()
    c=io.read_json(args.config)
    io.required(c,["table","delimiter","id_column","age_column","education_column","min_tests","session"])
    src=(args.config.parent / c["table"]).resolve()
    with src.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f,delimiter=c["delimiter"]))
    ids=[r[c["id_column"]].strip() for r in rows]
    keys=list(c["tests"])
    if not rows or "" in ids or len(ids)!=len(set(ids)):
        raise ValueError("Fournir une ligne par sujet pour la session indiquée")
    if not 1 <= c["min_tests"] <= len(keys) or not set(c["tests"].values()) <= {-1,1}:
        raise ValueError("Vérifier le minimum de tests et leur orientation")
    missing={"",*map(str,c["missing_codes"])}
    def value(r,k):
        return float(r[k]) if r[k].strip() not in missing else np.nan
    age=np.array([value(r,c["age_column"]) for r in rows])
    edu=np.array([value(r,c["education_column"]) for r in rows])
    x=np.array([[value(r,k) for k in keys] for r in rows])
    if not np.isfinite(age).all() or np.isinf(edu).any() or np.isinf(x).any():
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
    folder=io.OUTPUT/"cohorte_cognition"/io.token(args.run)
    with io.record(folder,"Score cognitif global et ACP descriptive",c,{"table":src}) as data:
        io.write_table(data/"scores.csv",["subject","session","age","education","cognition","tests_used"],
                       zip(ids,[c["session"]]*len(ids),age,edu,score,counts))
        io.write_table(data/"standardisation.csv",["test","sign","mean","std","n"],zip(keys,c["tests"].values(),mean,std,np.isfinite(x).sum(axis=0)))
        pcs=[f"PC{i+1}" for i in range(len(s))]
        io.write_table(data/"acp_scores.csv",["subject",*pcs],([i,*v] for i,v in zip(ids,u*s)))
        io.write_table(data/"acp_axes.csv",["test",*pcs],([k,*v] for k,v in zip(keys,vt.T)))
        io.write_table(data/"acp_variance.csv",["component","fraction"],zip(pcs,variance))
        def draw(fig,ax):
            ax.plot(np.arange(1,len(s)+1),np.cumsum(variance),"o-",color="#4c72b0")
            ax.set(xlabel="Nombre de composantes",ylabel="Fraction de variance cumulée",title="ACP cognitive descriptive")
        io.figure(folder,"acp",draw)
    print(folder)


if __name__=="__main__":
    main()

