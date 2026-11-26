#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source config/cluster.env
python_bin="${PYTHON_BIN:-python3}"
# Une chaîne séquentielle par processus ; adapter le parallélisme au niveau des jobs.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
if (( $# < 3 )); then
    echo "Usage: bash script/lancer_1000brains.sh ETAPE CONFIG RUN [options de pipeline.py]" >&2
    exit 2
fi
mode=$1
cfg=$2
run=$3
shift 3
if [[ $mode == pilote ]]; then
    for step in controle statique m1 cibles; do
        "$python_bin" script/pipeline.py "$step" --config "$cfg" --run "$run"
    done
elif [[ $mode == jumeau ]]; then
    for step in simuler mfm tmfm evaluer; do
        "$python_bin" script/pipeline.py "$step" --config "$cfg" --run "$run" "$@"
    done
else
    exec "$python_bin" script/pipeline.py "$mode" --config "$cfg" --run "$run" "$@"
fi

