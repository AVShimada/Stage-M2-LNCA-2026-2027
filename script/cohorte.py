"""Traitements après validation du pilote. Une session et une version par analyse de cohorte."""
import argparse
import csv
from pathlib import Path

import numpy as np

import commun as io
from metadonnees import session_id, subject_id
from pipeline import STAGES

BASE = ("controle","statique","m1","cibles")


def roster(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        entries = list(csv.DictReader(f))
    if not entries:
        raise ValueError("Le manifeste de cohorte est vide")
    records, used = [], set()
    for row in entries:
        p = (path.parent / row["config"]).resolve()
        c = io.read_json(p)
        io.required(c,["subject","session","version"])
        key = (subject_id(c["subject"]), session_id(c["session"]), c["version"])
        if key in used:
            raise ValueError(f"Enregistrement en double : {key}")
        used.add(key)
        records.append((c,io.token(row["run"]),row,p))
    return records


def optional_number(row,key):
    text = (row.get(key) or "").strip()
    if text.casefold() in ("", "nan"):
        return np.nan
    v = float(text)
    if not np.isfinite(v):
        raise ValueError(f"Valeur non finie dans {key}")
    return v


def index_table(path, columns):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        missing = {"subject", "session", *columns} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Colonnes absentes de {path.name} : {sorted(missing)}")
        indexed = {}
        for row in reader:
            key = (subject_id(row["subject"]), session_id(row["session"]))
            if key in indexed:
                raise ValueError(f"Enregistrement en double dans {path.name} : {key}")
            indexed[key] = row
        return indexed


def get_metrics(c,run,files,idx):
    paths = {}
    for key in ("statique","m1","cibles"):
        out, _ = io.completed(c,STAGES[key][0],run)
        paths[key] = out/"data"
        files[f"{idx}_{key}"] = out/"manifest.json"
    static = io.read_json(paths["statique"]/"metrics.json")
    speed = io.read_json(paths["m1"]/"metrics.json")
    targets = io.read_json(paths["cibles"]/"targets.json")
    return [static["sc_mean"],static["fc_z_mean"],speed["Vtyp_mid_emp_0"],*targets.values()],paths


def ols(x,y):
    a = np.column_stack((np.ones(len(x)),x))
    if len(y)<=a.shape[1] or np.linalg.matrix_rank(a)<a.shape[1]:
        raise ValueError("Effectif ou rang insuffisant pour le modèle")
    return np.linalg.lstsq(a,y,rcond=None)[0]


def predict(x,b):
    return b[0]+x@b[1:]


def r2(y,yhat):
    sst = ((y-y.mean())**2).sum()
    if sst==0:
        raise ValueError("Cible constante")
    return float(1-((y-yhat)**2).sum()/sst)


def associations(out,a,seed,nboot,nperm,subjects):
    y = a[:,2]
    x = a[:,[0,1,3,4,5]]  # age, education, SC, FC, Vtyp
    models = {"age":[0],"age_edu":[0,1],"age_SC":[0,2],"age_FC":[0,3],"age_dFC":[0,4],
              "age_imagerie":[0,2,3,4],"complet":[0,1,2,3,4]}
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(len(y)),5)
    rows, predictions = [], []
    for name,inds in models.items():
        xx = x[:,inds]
        fit = predict(xx,ols(xx,y))
        cv = np.empty_like(y)
        for f,test in enumerate(folds):
            train = np.ones(len(y),dtype=bool)
            train[test] = False
            cv[test] = predict(xx[test],ols(xx[train],y[train]))
            predictions.extend((subjects[i],name,f,float(y[i]),float(cv[i])) for i in test)
        rows.append([name,len(y),r2(y,fit),r2(y,cv)])
    io.write_table(out/"data/modeles.csv",["modele","n","R2_ajustement","R2_hors_pli"],rows)
    io.write_table(out/"data/predictions.csv",["subject","modele","pli","observe","predit"],predictions)
    def delta(xx,yy):
        return r2(yy,predict(xx,ols(xx,yy)))-r2(yy,predict(xx[:,:2],ols(xx[:,:2],yy)))
    observed = delta(x,y)
    base = predict(x[:,:2],ols(x[:,:2],y))
    residual = y-base
    perm = [delta(x,base+rng.permutation(residual)) for _ in range(nperm)]
    boot = []
    for _ in range(nboot):
        idx = rng.integers(len(y),size=len(y))
        try:
            boot.append(delta(x[idx],y[idx]))
        except ValueError:
            continue
    if not boot:
        raise ValueError("Aucun bootstrap estimable")
    io.write_json(out/"data/apport_imagerie.json",dict(delta_R2=observed,bootstrap_IC95=np.percentile(boot,[2.5,97.5]).tolist(),
                  bootstrap_valides=len(boot),p_permutation=(1+sum(v>=observed for v in perm))/(1+nperm)))
    np.savez(out/"data/reechantillonnages.npz",bootstrap=boot,permutation=perm)
    def draw(fig,ax):
        ix=np.arange(len(rows))
        ax.bar(ix-0.18,[r[2] for r in rows],0.36,label="Ajustement")
        ax.bar(ix+0.18,[r[3] for r in rows],0.36,label="Hors pli")
        ax.set_xticks(ix,[r[0] for r in rows],rotation=35,ha="right")
        ax.set(ylabel="R²",title="Cognition · cohorte commune")
        ax.legend()
    io.figure(out,"modeles",draw)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=["individu","assembler","groupes","modeles"])
    p.add_argument("--table",required=True,type=Path)
    p.add_argument("--scores",type=Path,help="scores.csv issu de cognition.py, joint sur sujet + session")
    p.add_argument("--demographie",type=Path,help="demographie.csv issu de metadonnees.py, joint sur sujet + session")
    p.add_argument("--reseaux",action="store_true",help="Calculer les graphes des matrices moyennes des groupes d’âge")
    p.add_argument("--mc",action="store_true",help="Produire la MC moyenne et une FCD individuelle par tertile d’âge")
    p.add_argument("--run")
    p.add_argument("--row",type=int,help="Numéro de ligne de données, à partir de 1")
    p.add_argument("--bootstrap",type=int,default=1000)
    p.add_argument("--permutations",type=int,default=5000)
    args=p.parse_args()
    records=roster(args.table)
    if args.stage=="individu":
        from pipeline import execute
        if args.row is None or not 1<=args.row<=len(records):
            p.error("--row doit désigner une ligne du manifeste")
        c,run,_,path=records[args.row-1]
        for stage in BASE:
            folder=io.run_dir(c,STAGES[stage][0],run)
            if folder.exists():
                _,m=io.completed(c,STAGES[stage][0],run)
                for info in m["inputs"].values():
                    if io.digest(Path(info["path"])) != info["sha256"]:
                        raise ValueError(f"Entrée changée depuis {stage}")
                continue
            execute(argparse.Namespace(stage=stage,config=path,run=run))
        return
    versions={(session_id(c["session"]),c["version"]) for c,_,_,_ in records}
    if len(versions)!=1:
        raise ValueError("Sélectionner une session et une version pour chaque cohorte transversale")
    files={"table":args.table.resolve()}
    scores, demographics = {}, {}
    if args.scores:
        files["cognition"]=args.scores.resolve()
        scores = index_table(args.scores, ["age", "education", "cognition"])
    if args.demographie:
        files["demographie"]=args.demographie.resolve()
        demographics = index_table(args.demographie, ["age", "education", "sex"])
    numbers, ids, sexes, paths = [], [], [], []
    matches = {"demographie_apparies": 0, "cognition_apparies": 0,
               "age_discordant_entre_sources": 0, "education_discordante_entre_sources": 0}
    convention=None
    for i,(c,run,row,path) in enumerate(records):
        common=(c["sc_weights"],c["m1"]["speed_definition"],c["m1"]["max_volumes_for_speed"],c["tr_seconds"])
        if convention is not None and common!=convention:
            raise ValueError("Harmoniser TR, échelle SC et définition de Vtyp pour la cohorte")
        convention=common
        metrics,dirs=get_metrics(c,run,files,i)
        key = (subject_id(c["subject"]), session_id(c["session"]))
        demography = demographics.get(key, {})
        score = scores.get(key, {})
        matches["demographie_apparies"] += int(bool(demography))
        matches["cognition_apparies"] += int(bool(score))
        demo=[]
        for k in ("age","education","cognition"):
            v=optional_number(row,k)
            w = optional_number(score,k)
            d = optional_number(demography,k)
            if np.isfinite(d):
                if np.isfinite(w) and not np.isclose(d,w):
                    counter = "age_discordant_entre_sources" if k == "age" else "education_discordante_entre_sources"
                    matches[counter] += 1
                w = d
            if np.isfinite(w):
                if np.isfinite(v) and not np.isclose(v,w):
                    raise ValueError(f"Valeur {k} du manifeste discordante pour {key[0]}")
                v=w
            demo.append(v)
        numbers.append(demo+metrics)
        ids.append(key[0])
        sexes.append(demography.get("sex",""))
        paths.append(dirs)
    a=np.array(numbers)
    from dynamique import TARGETS
    cols=["age","education","cognition","SC_mean","FC_z_mean","Vtyp_mid",*TARGETS]
    out=io.stage_dir("cohorte_"+args.stage)/io.token(args.run)
    session, version = next(iter(versions))
    opts=dict(session=session,version=version,seed=42,
              bootstrap=args.bootstrap,permutations=args.permutations,reseaux=args.reseaux,mc=args.mc)
    with io.record(out,"Analyse de cohorte · "+args.stage,opts,files) as data:
        headers = ["subject", *cols]
        rows = ([i,*v] for i,v in zip(ids,a))
        if args.demographie:
            headers.append("sex")
            rows = ([*row, sex] for row, sex in zip(rows, sexes))
        io.write_table(data/"sujets.csv", headers, rows)
        if args.demographie or args.scores:
            io.write_json(data/"jointure_metadonnees.json", matches)
        if args.stage=="modeles":
            keep=np.all(np.isfinite(a[:,:6]),axis=1)
            if keep.sum()<10 or min(args.bootstrap,args.permutations)<1:
                raise ValueError("Effectif complet ou nombre de rééchantillonnages insuffisant")
            associations(out,a[keep],42,args.bootstrap,args.permutations,np.asarray(ids)[keep])
        elif args.stage=="groupes":
            groups={"moins55":a[:,0]<55,"55plus":a[:,0]>=55}
            labels0=None
            for c,run,_,_ in records:
                folder,_=io.completed(c,STAGES["controle"][0],run)
                labels=(folder/"data/roi_order.csv").read_bytes()
                if labels0 is not None and labels0!=labels:
                    raise ValueError("Harmoniser l'ordre des régions avant une moyenne de matrices")
                labels0=labels
            label_names=io.labels(io.source(records[0][0],"labels_bold"))
            grouped={}
            group_matrices = {}
            for name,sel in groups.items():
                if not sel.any():
                    continue
                matrices={}
                for key in ("sc","fc_fisher_z"):
                    average=sum(np.load(paths[i]["statique"]/(key+".npy")) for i in np.flatnonzero(sel))/sel.sum()
                    np.save(data/f"{key}_{name}.npy",average)
                    group_matrices.setdefault(key, {})[name] = average
                    matrices["SC" if key=="sc" else "FC"]=average
                if args.reseaux:
                    from graphes import graph_metrics
                    folder=out/name
                    for sub in ("data","figures"):
                        (folder/sub).mkdir(parents=True)
                    graph_metrics(folder,label_names,matrices,records[0][0]["graphs"],42)
                    grouped[name]={k:np.load(folder/f"data/{k}_consensus.npz")["labels"] for k in matrices}
            for key, means in group_matrices.items():
                limits = (min(m.min() for m in means.values()), max(m.max() for m in means.values()))
                for name, average in means.items():
                    io.heatmap(out, key+"_"+name, average, f"{key} · {name} (n={groups[name].sum()})", limits=limits)
            if len(grouped)==2:
                for key in ("FC","SC"):
                    young,old=grouped["moins55"][key],grouped["55plus"][key]
                    counts=np.zeros((int(young.max()),int(old.max())),dtype=int)
                    np.add.at(counts,(young-1,old-1),1)
                    io.write_table(data/f"transitions_{key}.csv",["roi","module_moins55","module_55plus"],zip(label_names,young,old))
                    io.heatmap(out,"transitions_"+key,counts,"Correspondance des modules · nombre de régions",
                               "Module (indice à partir de 0)",style="transitions_"+key.lower())
            if args.mc:
                valid=np.isfinite(a[:,0])
                if valid.sum()<3:
                    raise ValueError("Au moins trois âges renseignés pour les tertiles")
                q=np.quantile(a[valid,0],[1/3,2/3])
                tertile=np.searchsorted(q,a[:,0],side="right")
                io.write_json(data/"tertiles.json",{"bornes_age":q.tolist(),"egalites":"groupe supérieur"})
                members=[]
                for g in range(3):
                    idx=np.flatnonzero(valid & (tertile==g))
                    if not len(idx):
                        continue
                    mean=np.zeros((4950,4950),dtype=np.float64)
                    for i in idx:
                        mean+=np.load(paths[i]["m1"]/"mc_20s.npy",mmap_mode="r")
                    mean/=len(idx)
                    np.save(data/f"mc_tertile{g+1}.npy",mean.astype("float32"))
                    i=idx[np.argmin(np.abs(a[idx,0]-np.median(a[idx,0])))]
                    for j in idx:
                        members.append((ids[j],a[j,0],g+1,int(j==i)))
                    fcd=np.load(paths[i]["m1"]/"fcd_20s.npy")
                    io.heatmap(out,f"fcd_tertile{g+1}",fcd,f"FCD individuelle · {ids[i]} · {a[i,0]:g} ans","Fenêtre retenue",style="fcd_m1")
                    io.heatmap(out,f"mc_tertile{g+1}",mean[::5,::5],f"MC moyenne · tertile {g+1} (n={len(idx)})","Indice de lien (× 5)",style="mc_m1")
                io.write_table(data/"sujets_tertiles.csv",["subject","age","tertile","exemple_FCD"],members)
            for j,label in enumerate(cols[1:6],1):
                def draw(fig,ax):
                    points=ax.scatter(a[:,0],a[:,j],c=a[:,0],cmap="viridis",s=10,alpha=0.75)
                    fig.colorbar(points,ax=ax,label="Âge (ans)")
                    ax.set(xlabel="Âge (ans)",ylabel=label,title=label+" selon l’âge")
                io.figure(out,label+"_age",draw)
    print(out)


if __name__=="__main__":
    main()

