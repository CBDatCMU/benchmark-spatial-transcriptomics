# Xenium Benchmarking — CMU Lane Center HPC

Scripts for running the Xenium benchmarking pipeline on the cluster.

## Paper

Sergio Marco Salas et al. *Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows.* Nature Methods 22, 813–823 (2025).
https://doi.org/10.1038/s41592-025-02617-2

Original pipeline code: https://github.com/Moldia/Xenium_benchmarking

## Dataset

Human spinal cord (inactive) Xenium dataset from Zenodo:
https://doi.org/10.5281/zenodo.11120922

To download:

```
cd ~/benchmark-spatial-transcriptomics
bash download_data.sh
```

This will place the data at `data/example_spinal_chord_inactive/`. If that folder already exists, the download is skipped.

If automatic download fails, download `example_spinal_chord_inactive.zip` manually from the link above and extract it into `data/`.

## Setup
Clone this repository along with the origianl Xenium benchmarking pipeline:
 
```
git clone --recurse-submodules https://github.com/CBDatCMU/benchmark-spatial-transcriptomics.git
cd benchmark-spatial-transcriptomics

# Update the submodule
git submodule update --init --recursive
```

Install the conda environment:
 
```
# Make sure that the current shell session is able to use conda
# This must be run once per login session.
source ~/miniconda3/etc/profile.d/conda.sh

# create and activate the conda environment
conda env create --file Xenium_benchmarking/xenium_benchmarking.yml --prefix ./envs/xb
conda activate ./envs/xb
pip install -e Xenium_benchmarking/
```

## Usage

```
# Submit by SLURM
sbatch submit_xenium.sh

# Run directly
cd ~/benchmark-spatial-transcriptomics
conda activate ~/xenium_benchmark/envs/xb
python run_xenium_benchmark.py
```

Steps 2 (Baysor) and 5 (SpaGE) are skipped. Baysor requires Docker and SpaGE requires a scRNA-seq reference dataset.
