"""Point d'entrée cluster. Exécuter depuis la racine du projet."""
import argparse
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from commun import (EDGE, OUTPUT, array_file, completed, digest, figure, heatmap, load_inputs, read_json, record,
                    required, run_dir, source, token, write_json, write_table)
import dynamique as dyn

STAGES = {
    "accord-fc": ("00_accord_fc", "Comparaison par lien de la FC fournie et de la FC recalculée sur toutes les lignes du BOLD livré."),
    "controle": ("00_controle", "Contrôle du BOLD, de la SC, des régions et du calendrier des volumes."),
    "statique": ("01_sc_fc", "SC et FC individuelles, force et distribution des poids."),
    "m1": ("02_dynamique_m1", "FCD et MC à 20 s, vitesses multi-fenêtres et substituts du M1."),
    "cibles": ("03_cibles_2026", "Dix cibles d’ajustement et descripteurs indépendants de validation."),
    "graphes": ("04_graphes", "Communautés, hubs, participation, efficacité et clustering de la FC et de la SC."),
    "distance": ("05_distance", "FC moyenne selon la distance entre les centres des régions."),
    "rich-club": ("06_rich_club", "Rich club pondéré et normalisé par des graphes nuls."),
    "modules-mc": ("07_modules_mc", "Modules du graphe des liens à partir de la MC individuelle."),
    "comparer": ("08_v4_moins_v1", "Effet de version pour le même sujet et le même enregistrement."),
    "simuler": ("09_simulation", "Simulation de contrôle du MFM et extraction des dix cibles."),
    "mfm": ("10_ajuster_mfm", "Ajustement du MFM autonome au sujet."),
    "tmfm": ("11_ajuster_tmfm", "Ajustement du tMFM à couplage fluctuant."),
    "evaluer": ("12_evaluer_jumeaux", "Évaluation des modèles sur de nouvelles graines et sur FC, couplage et MC."),
}


def prior(c, name, run, files):
    path, manifest = completed(c, STAGES[name][0], run)
    for key, p in files.items():
        if key in manifest["inputs"] and digest(p) != manifest["inputs"][key]["sha256"]:
            raise ValueError(f"Entrée modifiée depuis {name} : {p}")
    return path, manifest


def inspect_inputs(cfg):
    print(f"{cfg['subject']} / {cfg['session']} / {cfg['version']} / {cfg['atlas']}")
    bold = array_file(source(cfg, "bold"))
    sc = array_file(source(cfg, "sc"))
    ref = array_file(source(cfg, "fc_reference_z"))
    print(f"Dimensions sur disque : BOLD {bold.shape}, SC {sc.shape}, FC Fisher z {ref.shape}")
    if cfg.get("ts_dimensions"):
        with source(cfg, "ts_dimensions").open(encoding="utf-8-sig", errors="replace") as f:
            print("TSdim (début du fichier) :", f.read(300).strip())
    if sc.shape != (100, 100) or not np.isfinite(sc).all() or not np.allclose(sc, sc.T, atol=1e-8):
        raise ValueError("SC à examiner : matrice 100 × 100 finie et symétrique attendue")
    if ref.shape != (100, 100) or not np.isfinite(ref[EDGE]).all() or not np.allclose(ref[EDGE], ref.T[EDGE], atol=1e-8):
        raise ValueError("FC à examiner : matrice 100 × 100, finie et symétrique hors diagonale attendue")
    if bold.ndim != 2 or 100 not in bold.shape:
        raise ValueError("La série BOLD doit posséder un axe de 100 régions")
    if bold.shape == (100, 100):
        orientation = cfg.get("bold_orientation")
        if orientation not in ("time_regions", "regions_time"):
            raise ValueError("Deux axes de 100 : renseigner l'orientation d'après la documentation")
    else:
        orientation = "time_regions" if bold.shape[1] == 100 else "regions_time"
        if cfg.get("bold_orientation") and cfg["bold_orientation"] != orientation:
            raise ValueError("L'orientation configurée contredit les dimensions de la série")
    x = bold if orientation == "time_regions" else bold.T
    print(f"bold_orientation = {orientation} ; {len(x)} échantillons dans le fichier")
    if len(x) < 3 or not np.isfinite(x).all() or np.any(np.ptp(x, axis=0) == 0):
        raise ValueError("Comparaison FC indisponible : vérifier censure, valeurs manquantes ou région constante")
    # Comparer sur l'échelle r pour que les diagonales et arctanh(1) n'interviennent pas.
    measured, supplied = dyn.fc(x)[EDGE], np.tanh(ref[EDGE])
    errors = np.abs(measured - supplied)
    print(f"FC Pearson sur toutes les lignes du fichier vs tanh(FC fournie), {len(errors)} liens")
    print(f"Écart absolu médian = {np.median(errors):.6g} ; maximum = {errors.max():.6g}")
    return measured, supplied


def fc_agreement(cfg, run):
    stage, title = STAGES["accord-fc"]
    out = run_dir(cfg,stage,run)
    inputs = {k: source(cfg,k) for k in ("bold","fc_reference_z","sc")}
    if cfg.get("ts_dimensions"):
        inputs["ts_dimensions"] = source(cfg,"ts_dimensions")
    with record(out,title,cfg,inputs) as data:
        measured, supplied = inspect_inputs(cfg)
        if np.ptp(supplied) == 0 or np.ptp(measured) == 0:
            raise ValueError("Un profil de liens est constant ; la corrélation entre profils est indéfinie")
        slope, intercept = np.polyfit(supplied,measured,1)
        delta = measured - supplied
        residual = measured - (slope * supplied + intercept)
        metrics = dict(n_links=len(delta), r_profiles=float(np.corrcoef(supplied,measured)[0,1]),
                       n_recalculated_higher=int(np.count_nonzero(delta>0)),
                       median_delta=float(np.median(delta)), max_abs_delta=float(np.max(np.abs(delta))),
                       slope=float(slope), intercept=float(intercept),
                       max_abs_residual=float(np.max(np.abs(residual))))
        write_json(data/"accord_fc.json",metrics)
        write_table(data/"liens.csv",["roi_i_0","roi_j_0","r_fourni","r_recalcule","difference","residu"],
                    zip(*EDGE,supplied,measured,delta,residual))
        def draw(fig, ax):
            axis = np.array([supplied.min(),supplied.max()])
            ax.axhline(0,color="0.35",linestyle="--",label="Égalité des deux FC")
            ax.scatter(supplied,delta,s=7,alpha=0.35,rasterized=True,color="#4c72b0")
            ax.plot(axis,(slope-1)*axis+intercept,color="red",label="Ajustement affine")
            ax.set(xlabel="FC fournie (r)",ylabel="FC recalculée − FC fournie (r)",
                   title=f"Écart par lien · {cfg['subject']} · {cfg['session']} · {cfg['version']}")
            ax.legend()
        figure(out,"accord_fc",draw)
    print(f"FC recalculée = {slope:.10g} × FC fournie {intercept:+.10g}")
    print(f"Écart résiduel absolu maximal = {metrics['max_abs_residual']:.6g}")
    print(out/"figures/accord_fc.png")


def execute(args):
    cfg = read_json(args.config)
    if args.stage == "entrees":
        inspect_inputs(cfg)
        return
    if args.stage == "champs":
        for k,v in cfg.items():
            if v in ("",None):
                print(k)
        return
    run = token(args.run)
    if args.stage == "accord-fc":
        fc_agreement(cfg,run)
        return
    if args.stage == "figures":
        gather(cfg,run)
        return
    stage, title = STAGES[args.stage]
    name = args.stage
    x, grid, n, sc, names, files = load_inputs(cfg)
    extra = dict(n_volumes_original=n, n_volumes_kept=len(grid), tr_seconds=cfg["tr_seconds"])
    tr = cfg["tr_seconds"]
    deps = {}
    prerequisites = {"statique":["controle"],"m1":["controle"],"cibles":["statique"],
                     "graphes":["statique"],"distance":["statique"],"rich-club":["statique"],
                     "modules-mc":["m1"],"comparer":["statique","m1","cibles"],"simuler":["cibles"],
                     "mfm":["cibles","simuler"],"tmfm":["cibles","simuler"],"evaluer":["cibles","mfm","tmfm"]}
    for dep in prerequisites.get(name, []):
        path, _ = prior(cfg,dep,run,files)
        deps[dep] = path / "data"
        files["source_"+dep] = path / "manifest.json"
    if name == "distance":
        files["coordinates_csv"] = source(cfg,"coordinates_csv")
    if name == "comparer":
        if not args.reference or not args.reference_run:
            raise ValueError("Fournir --reference et --reference-run pour la V1")
        ref = read_json(args.reference)
        if (ref["subject"],ref["session"],ref["version"],cfg["version"]) != (cfg["subject"],cfg["session"],"V1","V4"):
            raise ValueError("Comparaison attendue : un sujet, une session, V4 moins V1")
        xr, gr, nr, sr, labels_ref, fref = load_inputs(ref)
        if (tr,n,names,cfg["m1"],cfg["seed"]) != (ref["tr_seconds"],nr,labels_ref,ref["m1"],ref["seed"]):
            raise ValueError("Comparer sur le même calendrier, ordre des régions et profil M1")
        if not np.array_equal(grid,gr):
            raise ValueError("Masques différents : préparer deux configurations avec les mêmes volumes communs")
        if cfg["sc_weights"] != ref["sc_weights"]:
            raise ValueError("Les conventions de poids SC doivent correspondre")
        for dep in prerequisites["comparer"]:
            path, _ = prior(ref,dep,token(args.reference_run),fref)
            deps["ref_"+dep] = path / "data"
            files["V1_"+dep] = path / "manifest.json"
    if name in ("simuler","mfm","tmfm","evaluer"):
        import jumeau as twin
        if not args.model:
            raise ValueError("Fournir --model config/jumeau.json")
        model = read_json(args.model)
        scm, fm = twin.setup(cfg,model,names)
        for dep in prerequisites[name]:
            if dep in ("simuler","mfm","tmfm"):
                m = read_json(files["source_"+dep])
                old = m["model_configuration"]
                if dep == "simuler":
                    # The control run precedes choosing the search ranges and GA settings.
                    unused = {"bounds","population","generations","optimizer_seeds","validation_seeds"}
                else:
                    unused = {"validation_seeds"}
                if any(old.get(k) != v for k,v in model.items() if k not in unused) or any(m["inputs"][k]["sha256"] != digest(p) for k,p in fm.items()):
                    raise ValueError(f"Configuration du modèle ou SC modifiée depuis {dep}")
        extra["model_configuration"] = model
        files.update(fm)
        empirical = read_json(deps["cibles"] / "targets.json")
        target = np.array([empirical[k] for k in dyn.TARGETS])
    out = run_dir(cfg,stage,run)
    with record(out,title,cfg,files,extra) as data:
        if name == "controle":
            write_table(data / "retained_volumes.csv", ["volume_0","time_s"], zip(grid,grid*tr))
            write_table(data / "roi_order.csv", ["position_0","roi"],enumerate(names))
            write_json(data / "quality.json", dict(n_original=n,n_kept=len(grid),continuous_blocks=[len(b) for b in dyn.blocks(grid)],
                       sc_min=float(sc.min()),sc_max=float(sc.max()),sc_weights=cfg["sc_weights"],tr_seconds=tr))
            def draw(fig,ax):
                full = np.full((n,3), np.nan)
                full[grid] = x[:,:3]
                ax.plot(np.arange(n)*tr,full,linewidth=0.8)
                ax.set(xlabel="Temps (s)",ylabel="Signal BOLD (unité fournie)",title="Trois premières régions · volumes retenus")
            figure(out,"controle_bold",draw)
        elif name == "statique":
            dyn.static(out,x,sc,names)
        elif name == "m1":
            dyn.m1(out,x,grid,tr,cfg)
        elif name == "cibles":
            dyn.empirical_targets(out,x,grid,tr,sc)
        elif name in ("graphes","distance","rich-club","modules-mc"):
            import graphes as graphs
            if name == "modules-mc":
                graphs.modules_mc(out,np.load(deps["m1"]/"mc_20s.npy",mmap_mode="r"),cfg["seed"],names)
            else:
                z = np.load(deps["statique"]/"fc_fisher_z.npy")
                if name == "graphes":
                    graphs.graph_metrics(out,names,{"FC":z,"SC":sc},cfg["graphs"],cfg["seed"])
                elif name == "distance":
                    graphs.distance(out,files["coordinates_csv"],names,z)
                else:
                    graphs.rich(out,z,cfg["graphs"]["rich_club_nulls"],cfg["seed"])
        elif name == "comparer":
            rows = []
            for dep, fname in (("statique","metrics.json"),("m1","metrics.json"),("cibles","targets.json"),("cibles","validation.json")):
                a,b = (read_json(deps[p+dep]/fname) for p in ("ref_",""))
                rows.extend([key,a[key],value,value-a[key]] for key,value in b.items())
            write_table(data/"differences.csv",["metric","V1","V4","V4_minus_V1"],rows)
            for key in ("fc_r","sc"):
                delta = np.load(deps["statique"]/(key+".npy")) - np.load(deps["ref_statique"]/(key+".npy"))
                np.save(data/(key+"_difference.npy"),delta)
                bound = float(np.abs(delta).max()) or 1.0
                heatmap(out,key+"_difference",delta,key+" · V4 moins V1",limits=(-bound, bound))
        elif name == "simuler":
            required(model["control_parameters"],["G","sigma"])
            y = twin.simulate(scm,cfg,model,model["control_parameters"],n,model["fit_seeds"][0])
            np.save(data/"simulated_bold.npy",y)
            v, _, _, fcd, *_ = dyn.feature_set(y[grid],grid,tr)
            write_json(data/"targets.json",dict(zip(dyn.TARGETS,v.tolist())))
            heatmap(out,"fcd_simulee",fcd,"FCD simulée · paramètres de contrôle","Fenêtre retenue",style="fcd_article")
        elif name in ("mfm","tmfm"):
            twin.fit(out,name,cfg,model,scm,target,grid,n)
        elif name == "evaluer":
            fits = {k:deps[k]/"best.json" for k in ("mfm","tmfm")}
            empirical.update(read_json(deps["cibles"]/"validation.json"))
            twin.evaluate(out,cfg,model,scm,sc,fits,empirical,grid,n)
    print(out)


def gather(c, run):
    out = OUTPUT/"exports"
    out.mkdir(parents=True,exist_ok=True)
    at = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    dest = out/(run_dir(c,"",run).name+"_figures_"+at+".zip")
    sources = []
    for name,(stage,_) in STAGES.items():
        folder = run_dir(c,stage,run)
        if folder.exists():
            state = folder / "manifest.json"
            if not state.is_file() or read_json(state).get("status") != "complete":
                print(f"Étape en attente ou incomplète, omise de l'export : {stage}",file=sys.stderr)
                continue
            _, m = completed(c,stage,run)
            sources.extend((stage,p) for p in [folder/"README.md",folder/"manifest.json",*sorted((folder/"figures").glob("*"))])
    if not sources:
        raise ValueError("Aucune étape terminée pour cet identifiant")
    with zipfile.ZipFile(dest,"x",compression=zipfile.ZIP_DEFLATED) as z:
        for stage,p in sources:
            z.write(p,str(Path(stage)/p.relative_to(run_dir(c,stage,run))))
    print(dest)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=["champs","entrees",*STAGES,"figures"])
    p.add_argument("--config",required=True)
    p.add_argument("--run",help="Identifiant unique, commun aux étapes de cette exécution")
    p.add_argument("--model")
    p.add_argument("--reference",help="Configuration V1")
    p.add_argument("--reference-run")
    execute(p.parse_args())


if __name__ == "__main__":
    try:
        main()
    except (ValueError,OSError,KeyError) as exc:
        sys.exit(str(exc))

