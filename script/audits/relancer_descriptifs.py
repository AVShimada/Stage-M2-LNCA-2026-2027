"""Relancer les descriptifs de la cohorte filtrée et ranger les anciens résultats."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys

import commun as io


def completed_runs(base):
    for path in base.glob("*/manifest.json"):
        state = io.read_json(path)
        if state["status"] == "complete":
            yield path.parent, state


def cohort_run(expected):
    candidates = [(state["started"], folder, state) for folder, state in
                  completed_runs(io.stage_dir("cohorte_tests_complets"))
                  if state["configuration"].get("filtre_tests_attendus") and state.get("n_sujets_retenus") == expected]
    if not candidates:
        raise ValueError(f"Aucun export terminé avec les tests attendus et {expected} sujets.")
    _, folder, state = max(candidates)
    return folder, state


def launch(script, *arguments):
    command = [sys.executable, str(io.ROOT / "script" / script), *arguments]
    subprocess.run(command, cwd=io.ROOT, check=True)


def cognitive_scores(run):
    config_path = io.ROOT / "config/cognition.json"
    cfg = io.read_json(config_path)
    source = (config_path.parent / cfg["table"]).resolve()
    source_hash = io.digest(source)
    code_hash = io.digest(io.ROOT / "script/cognition.py")
    candidates = [(state["started"], folder.name) for folder, state in completed_runs(io.stage_dir("cohorte_cognition"))
                  if state["configuration"] == cfg
                  and state.get("inputs", {}).get("table", {}).get("sha256") == source_hash
                  and state.get("code_sha256", {}).get("cognition.py") == code_hash
                  and (folder / "data/scores.csv").is_file()]
    if candidates:
        selected = max(candidates)[1]
        print(f"Scores cognitifs de référence réutilisés ({selected}).", flush=True)
        return selected
    launch("cognition.py", "--config", str(config_path), "--run", run)
    return run


def archive_previous(base, cohort_state):
    source_run = io.token(cohort_state["configuration"]["source_run"])
    origin = (io.stage_dir("cohorte_second_scan") / source_run / "data/second_scan_demographie.xlsx").resolve()
    previous_n = cohort_state["n_sujets_source"]
    archive = (io.ROOT / "archives" / f"descriptifs_N{previous_n}").resolve()
    if not archive.is_relative_to(io.ROOT.resolve()):
        raise ValueError("Le dossier d’archives doit se trouver dans le projet.")
    previous = []
    for folder, state in completed_runs(base):
        demo_path = state.get("inputs", {}).get("demographie", {}).get("path")
        if state.get("n_sujets") == previous_n and demo_path and Path(demo_path).resolve() == origin:
            previous.append(folder)
    for folder in previous:
        if folder.is_symlink() or folder.resolve().parent != base.resolve():
            raise ValueError("Dossier de résultats inattendu pour l’archivage.")
        destination = archive / folder.name
        if destination.exists() or destination.is_symlink():
            print(f"Archive déjà présente pour {folder.name}, résultat conservé à son emplacement actuel.")
            continue
        archive.mkdir(parents=True, exist_ok=True)
        folder.rename(destination)
        print(f"Archivé : {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effectif", type=int, default=409, help="Nombre de sujets dans l’export filtré attendu")
    args = parser.parse_args()
    cohort, state = cohort_run(args.effectif)
    run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + f"_N{args.effectif}"
    scores = cognitive_scores(run)
    launch("decrire_second_scan.py", "--cohorte-run", cohort.name, "--scores-run", scores, "--run", run)

    base = io.stage_dir("cohorte_second_scan_descriptif")
    result = base / run
    summary = io.read_json(result / "manifest.json")
    if (summary["status"] != "complete" or summary["n_sujets"] != args.effectif
            or summary["n_scores_cognitifs"] != args.effectif):
        raise ValueError("L’effectif démographique ou cognitif du descriptif diffère de l’effectif attendu.")
    archive_previous(base, state)
    print(f"Descriptifs vérifiés pour {args.effectif} sujets : {result}")


if __name__ == "__main__":
    main()
