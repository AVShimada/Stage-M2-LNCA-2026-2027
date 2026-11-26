#!/usr/bin/env bash
# Reprendre une relance interrompue à partir du rich club, sur le nœud initial.
set -euo pipefail
cd "$(dirname "$0")/../.."
run=${1:?Donner l’identifiant UTC du lancement}
[[ $run =~ ^[0-9]{8}T[0-9]{6,12}Z$ ]] || { echo "Identifiant UTC invalide" >&2; exit 2; }
tag=${run}_sub-0001_ses-1_V4
export PYTHON_BIN="${PYTHON_BIN:-$HOME/.venvs/lnca-1000brains/bin/python}"
export PYTHONPATH="$PWD/script" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export_dir="output/04_exports/$run"
cfg="$export_dir/pilote.json"
[[ -f "$cfg" ]] || { echo "Fiche du lancement introuvable : $cfg" >&2; exit 1; }
exec 9>"$export_dir/reprise.lock"
flock -n 9 || { echo "Une reprise est déjà active."; exit 1; }
if pgrep -au "$(id -u)" -f "[p]ython[^ ]* .*relancer_workflow[.]py|[p]ython[^ ]* .*script/.*[.]py .*$run"; then
    echo "Un lanceur ou une étape est encore actif sur $(hostname). Vérifier le processus affiché."
    exit 1
fi
rich="output/02_pilote/10_rich_club/$tag"
field() { "$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2]) or "")' "$1" "$2"; }
first_manifest="output/02_pilote/09_distance/$tag/manifest.json"
job=$(field "$first_manifest" job_id)
if [[ -f "$rich/manifest.json" ]]; then
    rich_job=$(field "$rich/manifest.json" job_id)
    job=${rich_job:-$job}
fi
if [[ -n "$job" ]]; then
    echo "Vérification du job d’origine $job"
    jobs=$(squeue -r -h -u "$(id -un)" -o '%i %F')
    if awk -v job="$job" '$1==job || $2==job {found=1} END {exit !found}' <<<"$jobs"; then
        echo "Le job $job est encore dans la file. Vérifier son état sur son nœud avant de reprendre."
        exit 1
    fi
fi
exec > >(tee -a "$export_dir/reprise.log") 2>&1
"$PYTHON_BIN" - "$cfg" "$run" <<'PY'
import sys
from commun import completed, read_json
cfg = read_json(sys.argv[1])
for key in ("01_sc_fc", "02_dynamique_m1", "03_cibles_2026", "05_distance"):
    completed(cfg, key, sys.argv[2])
PY
stage() {
    local folder=$1 status saved
    shift
    if [[ -e "$folder" ]]; then
        status=$(field "$folder/manifest.json" status)
        case "$status" in
            complete) echo "Déjà terminé : $folder"; return ;;
            running|failed)
                mkdir -p "archives/reprises/$tag"
                saved=$(mktemp -d "archives/reprises/$tag/$(basename "$(dirname "$folder")").XXXXXX")
                mv -T -- "$folder" "$saved/$tag"
                echo "Tentative précédente conservée dans $saved/$tag" ;;
            *) echo "Statut inattendu dans $folder : $status"; exit 1 ;;
        esac
    fi
    echo "Calcul : $folder"
    "$@"
}
pipeline() { "$PYTHON_BIN" script/pipeline.py "$1" --config "$cfg" --run "$run"; }
stage "$rich" pipeline rich-club
stage "output/02_pilote/11_modules_mc/$tag" pipeline modules-mc
simulation="output/03_jumeau/01_simulation/$tag"
stage "$simulation" "$PYTHON_BIN" script/audits/controle_jumeau.py --source-run "$run" --run "$run" --simulation-seule
stage "output/03_jumeau/02_comparaison/$tag" "$PYTHON_BIN" script/audits/comparer_jumeau.py --source-run "$simulation" --run "$run"
pipeline figures
bash script/architecture_tree.sh
echo "Terminé. Archive des figures dans $export_dir"
