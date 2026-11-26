"""Rejouer les descriptifs, le pilote V4, ses graphes et la comparaison du jumeau."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shlex
import subprocess
import sys

import commun as io


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Fiche du pilote, par défaut celle validée le 20 novembre")
    args = parser.parse_args()
    config = args.config
    if config is None:
        config = io.ROOT / "config/pilote_ordre_M1_20261120T150810Z.json"
        if not config.exists():
            config = io.ROOT / "config/pilote.json"
    config = config.resolve()
    cfg = io.read_json(config)
    if (cfg["subject"], cfg["session"], cfg["version"]) != ("sub-0001", "ses-1", "V4"):
        parser.error("Cette relance concerne sub-0001 / ses-1 / V4.")
    # Ces références anatomiques sont celles du pilote déjà exécuté sur le cluster.
    for key in ("labels_bold", "labels_sc"):
        if not cfg.get(key):
            cfg[key] = str(io.ROOT / "atlas/labels_M1.txt")
            cfg["roi_order_note"] = "Pilote technique : ordre commun repris du M1, à confirmer."
    if not cfg.get("coordinates_csv"):
        cfg["coordinates_csv"] = str(io.ROOT / "atlas/coordonnees_M1.csv")

    io.load_inputs(cfg)
    for key in ("coordinates_csv", "fc_reference_z"):
        io.source(cfg, key)
    import graphes  # Vérifier les dépendances de graphes avant les calculs de population.
    import numba
    import openpyxl

    run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    pilot_tag = io.run_dir(cfg, "", run).name
    export = io.stage_dir("exports") / run
    export.mkdir(parents=True, exist_ok=False)
    active = export / "pilote.json"
    io.write_json(active, cfg)
    index = io.OUTPUT / "README.md"
    if not index.exists():
        index.write_text("# Sorties 1000BRAINS\n\nChaque analyse conserve son manifeste et ses entrées. "
                         "Un même identifiant UTC relie les étapes d’une relance.\n\n"
                         + "".join(f"- `{path}/`\n" for path in io.STAGE_PATHS.values()), encoding="utf-8")
    env = dict(os.environ, PYTHONPATH=str(io.ROOT / "script"), PYTHON_BIN=sys.executable,
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
    print(f"RUN {run}\nFiche pilote : {config}\nJournal : {export / 'relance.log'}", flush=True)

    with (export / "relance.log").open("x", encoding="utf-8") as log:
        def launch(script, *options):
            program = "bash" if script.endswith(".sh") else sys.executable
            command = [program, f"script/{script}", *map(str, options)]
            print(f"\n{shlex.join(command)}", flush=True)
            log.write(f"\n{shlex.join(command)}\n")
            log.flush()
            with subprocess.Popen(command, cwd=io.ROOT, env=env, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, encoding="utf-8", errors="replace") as proc:
                for line in proc.stdout:
                    print(line, end="", flush=True)
                    log.write(line)
                    log.flush()
                if proc.wait():
                    raise SystemExit(f"Échec de {script}. Journal : {log.name}")

        launch("metadonnees.py", "--run", run)
        launch("cognition.py", "--config", "config/cognition.json", "--run", run)
        launch("audits/manquants.py", "--run", run)
        launch("extraire_second_scan.py", "--run", run)
        launch("audits/filtrer_second_scan.py", "--source-run", run, "--run", run)
        launch("decrire_second_scan.py", "--population", "filtree", "--cohorte-run", run,
               "--scores-run", run, "--run", run)
        stats = io.read_json(io.stage_dir("cohorte_second_scan_descriptif") / run / "manifest.json")
        if (stats["n_sujets"], stats["n_scores_cognitifs"]) != (409, 409):
            raise SystemExit(f"Effectif observé : {stats['n_sujets']} sujets et {stats['n_scores_cognitifs']} "
                             f"scores dans la cohorte attendue de 409. Voir {log.name}")

        def pipeline(step):
            launch("pipeline.py", step, "--config", active, "--run", run)

        pipeline("controle")
        pipeline("accord-fc")
        launch("apercu_sc_fc.sh", active, run)
        pipeline("statique")
        pipeline("m1")
        launch("comparer_fcd_v1_v4.sh", active, run)
        pipeline("cibles")
        for stage in ("graphes", "distance", "rich-club", "modules-mc"):
            pipeline(stage)
        launch("audits/controle_jumeau.py", "--source-run", run, "--run", run, "--simulation-seule")
        control = io.stage_dir("09_jumeau_controle_technique") / pilot_tag
        launch("audits/comparer_jumeau.py", "--source-run", control, "--run", run)
        pipeline("figures")
        launch("architecture_tree.sh")
    print(f"Relance terminée. Comparaisons : {io.stage_dir('10_jumeau_comparaison') / pilot_tag / 'figures'}")
    print(f"Archive des figures : {export}")


if __name__ == "__main__":
    main()
