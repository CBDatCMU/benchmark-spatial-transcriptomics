#!/bin/bash
#SBATCH --job-name=xenium_benchmark
#SBATCH --output=logs/slurm_%x_%j.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --partition=short1
 
echo "Job started : $(date)"
echo "Node        : $(hostname)"
echo "Node list   : $SLURM_NODELIST"
echo "CPUs        : $SLURM_CPUS_PER_TASK"
echo "Memory      : $SLURM_MEM_PER_NODE MB"
 
# ── Activate conda env ────────────────────────────────────────────────────────
export PATH=/shared/containers/images/miniconda3/bin:$PATH
source /shared/containers/images/miniconda3/etc/profile.d/conda.sh
conda activate ~/benchmark-spatial-transcriptomics/envs/xb
 
# ── Run pipeline ──────────────────────────────────────────────────────────────
cd ~/benchmark-spatial-transcriptomics/Xenium_benchmarking
python ~/benchmark-spatial-transcriptomics/run_xenium_benchmark.py
 
echo "Job finished: $(date)"
