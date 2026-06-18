import os
import matplotlib
matplotlib.use("Agg")

import scanpy as sc
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import matplotlib.pyplot as plt

import cell2location
import torch

if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    
    visium_dir = "C:/Users/Gabe/Desktop/Visium"
    
    results_folder = (
       "C:/Users/Gabe/Desktop/multiome_poa/Python/cell2location/results_2026_06_17"
    )
    
    ref_run_name = f"{results_folder}/reference_signatures"
    
    run_name = f"{results_folder}/cell2location_map"
    
    os.makedirs(run_name, exist_ok=True)
    
    # ------------------------------------------------------------------
    # Reconstruct AnnData from MTX + CSV
    # ------------------------------------------------------------------
    
    print("Loading sketch matrix...")
    
    X = mmread(
        f"{visium_dir}/sketch_counts.mtx"
    ).T.tocsr()
    
    obs = pd.read_csv(
        f"{visium_dir}/sketch_meta.csv",
        index_col=0
    )
    
    var = pd.read_csv(
        f"{visium_dir}/sketch_genes.csv"
    )
    
    var.index = var["gene"]
    
    adata_vis = sc.AnnData(
        X=X,
        obs=obs,
        var=var
    )
    
    print(adata_vis)
    
    adata_vis.write_h5ad(
        f"{visium_dir}/sketch_visium.h5ad"
    )
    
    print("Saved sketch_visium.h5ad")
    
    # ------------------------------------------------------------------
    # Load reference signatures
    # ------------------------------------------------------------------
    
    print("Loading reference signatures...")
    
    inf_aver = pd.read_csv(
        f"{ref_run_name}/reference_signatures.csv",
        index_col=0
    )
    
    print(
        "Reference signature matrix:",
        inf_aver.shape
    )
    
    print(
        "Reference cell types:",
        len(inf_aver.columns)
    )
    
    # ------------------------------------------------------------------
    # Intersect genes
    # ------------------------------------------------------------------
    
    shared_genes = adata_vis.var_names.intersection(
        inf_aver.index
    )
    
    print(
        f"Shared genes: "
        f"{len(shared_genes)} / {len(inf_aver)}"
    )
    
    adata_vis = adata_vis[:, shared_genes].copy()
    
    inf_aver = inf_aver.loc[shared_genes]
    
    # ------------------------------------------------------------------
    # Batch column
    # ------------------------------------------------------------------
    
    print("\nMetadata columns:")
    print(adata_vis.obs.columns.tolist())
    
    BATCH_KEY = "Slice"
    
    if BATCH_KEY not in adata_vis.obs.columns:
        raise KeyError(
            f"{BATCH_KEY} not found.\n"
            f"Available columns:\n"
            f"{adata_vis.obs.columns.tolist()}"
        )
    
    print(
        f"\nBatch key: {BATCH_KEY}"
    )
    
    print(
        adata_vis.obs[BATCH_KEY].value_counts()
    )
    
    # ------------------------------------------------------------------
    # Setup cell2location
    # ------------------------------------------------------------------
    
    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis,
        batch_key=BATCH_KEY
    )
    
    # ------------------------------------------------------------------
    # Hyperparameters
    # ------------------------------------------------------------------
    
    N_CELLS = 3
    
    DETECTION_ALPHA = 200
    
    mod_sp = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=N_CELLS,
        detection_alpha=DETECTION_ALPHA,
    )
    
    mod_sp.view_anndata_setup()
    
    print("\nCUDA status")
    
    print(
        "CUDA available:",
        torch.cuda.is_available()
    )
    
    if torch.cuda.is_available():
        print(
            "GPU:",
            torch.cuda.get_device_name(0)
        )
    
    # ------------------------------------------------------------------
    # Train
    # ------------------------------------------------------------------
    
    mod_sp.train(
        max_epochs=1000,
        batch_size=2048,
        train_size=1,
        accelerator="gpu",
        datasplitter_kwargs={
            "num_workers": 8,
            "persistent_workers": True,
            "pin_memory":True
        },
    )
    
    # ------------------------------------------------------------------
    # Export posterior
    # ------------------------------------------------------------------
    
    adata_vis = mod_sp.export_posterior(
        adata_vis,
        sample_kwargs={
            "num_samples": 100,
            "batch_size": 512,
            "accelerator": "cpu",
        },
    )
    
    # ------------------------------------------------------------------
    # Save model
    # ------------------------------------------------------------------
    
    mod_sp.save(
        run_name,
        overwrite=True
    )
    
    adata_vis.write_h5ad(
        f"{visium_dir}/adata_vis_sketch_trained.h5ad"
    )
    
    # ------------------------------------------------------------------
    # Save abundance estimates
    # ------------------------------------------------------------------
    
    adata_vis.obsm[
        "q05_cell_abundance_w_sf"
    ].to_csv(
        f"{run_name}/abundance_q05.csv"
    )
    
    adata_vis.obsm[
        "means_cell_abundance_w_sf"
    ].to_csv(
        f"{run_name}/abundance_means.csv"
    )
    
    adata_vis.obsm[
        "q95_cell_abundance_w_sf"
    ].to_csv(
        f"{run_name}/abundance_q95.csv"
    )
    
    # ------------------------------------------------------------------
    # QC plots
    # ------------------------------------------------------------------
    
    mod_sp.plot_history(1000)
    
    plt.savefig(
        f"{run_name}/sp_train_history.png",
        bbox_inches="tight"
    )
    
    plt.close("all")
    
    mod_sp.plot_QC()
    
    plt.savefig(
        f"{run_name}/sp_qc.png",
        bbox_inches="tight"
    )
    
    plt.close("all")
    
    # ------------------------------------------------------------------
    # Region summaries
    # ------------------------------------------------------------------
    
    cell_types = (
        adata_vis.uns["mod"]["factor_names"]
    )
    
    adata_vis.obs[
        cell_types
    ] = adata_vis.obsm[
        "means_cell_abundance_w_sf"
    ].values
    
    region_abundance = (
        adata_vis.obs[
            ["dissection"] + list(cell_types)
        ]
        .groupby("dissection")
        .mean()
    )
    
    region_abundance.to_csv(
        f"{run_name}/abundance_by_dissection_region.csv"
    )
    
    print("\nFinished.")
    print(region_abundance)
    print(f"\nOutputs saved to:\n{run_name}")
