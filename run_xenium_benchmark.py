#!/usr/bin/env python

import os
import sys

# ── Banksy submodule path ──────────────────────────────────────────────────────
_repo = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_repo, 'Xenium_benchmarking'))
sys.path.insert(0, os.path.join(_repo, 'Xenium_benchmarking', 'Banksy_py'))

import time
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for HPC (no display)
import matplotlib.pyplot as plt
import squidpy as sq

# ── Banksy submodule path ──────────────────────────────────────────────────────
banksy_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'Banksy_py')
)
if banksy_path not in sys.path:
    sys.path.insert(0, banksy_path)

# ── xb package imports ─────────────────────────────────────────────────────────
from xb.formatting import *
from xb.plotting import *
from xb.preprocessing import *
from xb.Spage_main import *
from xb.calculating import *
from xb.domain_identification import *
from xb.neighborhood import *
import xb.calculating
xb.calculating.sq = sq

# ── Logging setup ──────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger(__name__)

# ── Timing helper ──────────────────────────────────────────────────────────────
_timings = {}

def timed(label):
    """Context manager that logs wall-clock time for a pipeline step."""
    class _Timer:
        def __enter__(self):
            self._start = time.time()
            log.info(f"START  {label}")
            return self
        def __exit__(self, *args):
            elapsed = time.time() - self._start
            _timings[label] = elapsed
            log.info(f"END    {label}  ({elapsed:.1f}s / {elapsed/60:.2f}min)")
    return _Timer()

# ══════════════════════════════════════════════════════════════════════════════
# 0. Parameters
# ══════════════════════════════════════════════════════════════════════════════
BASE        = os.path.expanduser('~/xenium_benchmark')
DATA_DIR    = os.path.join(BASE, 'data')
REPO_DIR    = os.path.join(BASE, 'Xenium_benchmarking')
OUTPUT_PATH = os.path.join(BASE, 'pipeline_output') + '/'
PLOT_PATH   = os.path.join(OUTPUT_PATH, 'figures') + '/'
LOG_PATH    = os.path.join(BASE, 'logs')

os.makedirs(OUTPUT_PATH, exist_ok=True)
os.makedirs(PLOT_PATH,   exist_ok=True)
os.makedirs(LOG_PATH,    exist_ok=True)

# Add file handler to logging now that LOG_PATH exists
fh = logging.FileHandler(os.path.join(LOG_PATH, 'benchmark.log'))
fh.setFormatter(logging.Formatter('%(asctime)s  %(levelname)s  %(message)s',
                                   datefmt='%Y-%m-%d %H:%M:%S'))
log.addHandler(fh)

SAMPLE_NAME = 'example_spinal_chord_inactive'
files       = [os.path.join(DATA_DIR, SAMPLE_NAME)]

save = True

max_nucleus_distance = 10
min_quality          = 0

clustering_params = {
    'normalization_target_sum': 100,
    'min_counts_x_cell':        40,
    'min_genes_x_cell':         15,
    'scale':                    False,
    'clustering_alg':           'louvain',
    'resolutions':              [0.2, 0.5, 1.1],
    'n_neighbors':              15,
    'umap_min_dist':            0.1,
    'n_pcs':                    0,
}

banksy_params = {
    'resolutions':       [0.9],
    'pca_dims':          [20],
    'lambda_list':       [0.8],
    'k_geom':            15,
    'max_m':             1,
    'nbr_weight_decay':  'scaled_gaussian',
    'cluster_algorithm': 'leiden',
}

hyperparameters_nbd = {
    'key':                  'louvain_1.1',
    'resolution':           1.0,
    'neighbors':            15,
    'clustering_algorithm': 'leiden',
}

hyperparameters_rbd = {
    'key':                  'louvain_1.1',
    'resolution':           0.2,
    'neighbors':            15,
    'clustering_algorithm': 'leiden',
}

log.info("Pipeline parameters set.")
log.info(f"  Sample : {SAMPLE_NAME}")
log.info(f"  Output : {OUTPUT_PATH}")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Format Xenium output to AnnData
# ══════════════════════════════════════════════════════════════════════════════
with timed("Step 1: Format Xenium → AnnData"):
    adata_raw = format_xenium_adata_final(
        files[0],
        SAMPLE_NAME,
        OUTPUT_PATH,
        use_parquet=True,
        save=False,
    )

    # Fix bytes columns (may occur on some systems)
    spots = adata_raw.uns['spots'].copy()
    for col in spots.columns:
        if spots[col].dtype == object:
            try:
                spots[col] = spots[col].apply(
                    lambda x: x.decode('utf-8') if isinstance(x, bytes) else x
                )
            except Exception:
                pass
    adata_raw.uns['spots'] = spots

    adata = keep_nuclei_and_quality(
        adata_raw,
        tag=SAMPLE_NAME,
        max_nucleus_distance=max_nucleus_distance,
        min_quality=min_quality,
        save=False,
        output_path=OUTPUT_PATH,
    )

    adata.obs_names_make_unique()
    adata.obs['expressed_genes'] = np.sum(adata.X > 0, axis=1)
    adata.obs['n_counts']        = np.sum(adata.X,     axis=1)
    adata.obs['n_genes']         = np.sum(adata.X > 0, axis=1)

    log.info(f"  adata shape after filtering: {adata.shape}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Baysor — SKIPPED
# ══════════════════════════════════════════════════════════════════════════════
log.info("Step 2: Baysor — skipped.")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Preprocess and cluster
# ══════════════════════════════════════════════════════════════════════════════
with timed("Step 3: Preprocess + cluster"):
    adata = preprocess_adata(
        adata,
        save=save,
        clustering_params=clustering_params,
        output_path=OUTPUT_PATH,
    )

# ══════════════════════════════════════════════════════════════════════════════
# 4. Domain identification
# ══════════════════════════════════════════════════════════════════════════════

# ── 4.1 Banksy ────────────────────────────────────────────────────────────────
with timed("Step 4.1: Banksy domain identification"):
    import random
    random_seed = 1234
    np.random.seed(random_seed)
    random.seed(random_seed)

    adata = sc.read(OUTPUT_PATH + 'combined_processed.h5ad')
    adata, adata_banksy = domains_by_banksy(
        adata,
        plot_path=PLOT_PATH,
        banksy_params=banksy_params,
        save=save,
    )

# ── 4.2 NBD ───────────────────────────────────────────────────────────────────
with timed("Step 4.2: Neighbors-based domain identification (NBD)"):
    adata.obs_names_make_unique()
    adata.obs['unique_cell_id'] = (
        adata.obs['sample'].astype(str) + '_' +
        adata.obs['cell_id'].astype(str)
    )
    adata, adata_nbd = domains_by_nbd(
        adata,
        hyperparameters_nbd=hyperparameters_nbd,
    )

# ── 4.3 RBD ───────────────────────────────────────────────────────────────────
with timed("Step 4.3: Read-based domain identification (RBD)"):
    adata, adata_rbd = domains_by_rbd(
        adata,
        hyperparameters_rbd=hyperparameters_rbd,
    )

# ── 4.4 Compare ───────────────────────────────────────────────────────────────
with timed("Step 4.4: Compare domain algorithms"):
    domain_keys = ['rbd_domain', 'nbd_domain', 'banksy_domain']
    ARI = compare_domains(
        adata,
        domain_keys=domain_keys,
        plot_path=PLOT_PATH,
        save=save,
    )

# ══════════════════════════════════════════════════════════════════════════════
# 5. SpaGE imputation — SKIPPED (no scRNA-seq reference)
# ══════════════════════════════════════════════════════════════════════════════
log.info("Step 5: SpaGE imputation — skipped (no scRNA-seq reference).")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Spatially variable genes (Moran's I)
# ══════════════════════════════════════════════════════════════════════════════
with timed("Step 6: Spatially variable genes (Moran's I)"):
    import squidpy as sq
    adata, hs_results = svf_moranI(adata, radius=50.0)
    log.info(f"  Top 10 SVGs:\n{hs_results.head(10)}")
    if save:
        hs_results.to_csv(OUTPUT_PATH + 'svf_moranI_results.csv')

# ══════════════════════════════════════════════════════════════════════════════
# 7. Neighborhood enrichment
# ══════════════════════════════════════════════════════════════════════════════
with timed("Step 7: Neighborhood enrichment (Squidpy)"):
    adata = nhood_squidpy(
        adata,
        sample_key='sample',
        radius=50.0,
        cluster_key='louvain_1.1',
        plot_path=PLOT_PATH,
        cmap='Blues',
    )

# ══════════════════════════════════════════════════════════════════════════════
# Timing summary
# ══════════════════════════════════════════════════════════════════════════════
log.info("\n" + "="*60)
log.info("BENCHMARK TIMING SUMMARY")
log.info("="*60)
total = 0
for step, secs in _timings.items():
    log.info(f"  {step:<45}  {secs/60:6.2f} min")
    total += secs
log.info("-"*60)
log.info(f"  {'TOTAL':<45}  {total/60:6.2f} min")
log.info("="*60)
