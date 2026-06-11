#!/bin/bash
#SBATCH --job-name=xenium_benchmark
#SBATCH --output=/home/bonniela/xenium_benchmark/logs/slurm_%j.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --partition=short1

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"

# ── Activate conda env ────────────────────────────────────────────────────────
export PATH=/containers/images/miniconda3/bin:$PATH
source /containers/images/miniconda3/etc/profile.d/conda.sh
conda activate ~/xenium_benchmark/envs/xb

# ── Run pipeline ──────────────────────────────────────────────────────────────
cd ~/xenium_benchmark/Xenium_benchmarking
python ~/xenium_benchmark/run_xenium_benchmark.py

echo "Job finished: $(date)"
