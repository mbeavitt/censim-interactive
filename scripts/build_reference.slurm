#!/bin/bash
# Build the real-data reference on CSD3. Array job: one task per (accession,chr)
# HOR table -> a partial .npz; merge afterwards into reference.dat.
#
# Setup + submit (from the dir containing build_reference.py):
#   mkdir -p logs partials
#   find <TRASH_ROOT> -name 'HORs_method_1_aligned_*_Chr*_CEN178.fasta_t_7_c_3.csv' \
#       | sort > pairs.txt
#   sbatch --array=1-$(wc -l < pairs.txt)%64 build_reference.slurm
#   # when the array finishes:
#   python build_reference.py merge --inputs partials --out reference.dat
#   # then scp reference.dat back next to the censim binary.
#
#SBATCH -J censim_ref
#SBATCH -o logs/ref_%A_%a.out
#SBATCH -e logs/ref_%A_%a.err
#SBATCH -t 00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -A CHANGE_ME_ACCOUNT          # your CSD3 project (e.g. XXX-SL3-CPU)
#SBATCH -p icelake                    # CPU work: use a CPU partition. (Wilkes3 =
                                      # ampere is GPU; no GPU is needed here.)

set -euo pipefail

# Adjust to your environment (module name varies on CSD3, or activate a venv):
module load python/3.11.0-icl 2>/dev/null || true
# python -m venv ~/censim_ref_env && source ~/censim_ref_env/bin/activate && pip install numpy

HOR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pairs.txt)
[ -n "$HOR" ] || { echo "no pair on line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$SCRIPT_DIR/build_reference.py" one \
    --hor "$HOR" \
    --out "partials/$(printf '%04d' "${SLURM_ARRAY_TASK_ID}").npz"
