import os
import matplotlib
matplotlib.use('Agg')  # must be first, before pyplot

import scanpy as sc
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import cell2location
from cell2location.models import RegressionModel
import torch


results_folder =    "C:/Users/Gabe/Desktop/Visium/POA"
ref_run_name   = f"{results_folder}/reference_signatures"
os.makedirs(ref_run_name, exist_ok=True)

# ── Load ───────────────────────────────────────────────────────────────
adata_ref = sc.read_h5ad("C:/Users/Gabe/Desktop/POA_ref.h5ad")

if "counts" in adata_ref.layers:
    adata_ref.X = adata_ref.layers["counts"].copy()
if sp.issparse(adata_ref.X):
    adata_ref.X = adata_ref.X.toarray()

sc.pp.filter_genes(adata_ref, min_cells=10)

# ── Setup ──────────────────────────────────────────────────────────────
RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key="pool",
    labels_key="sub.cluster",
)

# ── Train ──────────────────────────────────────────────────────────────
mod = RegressionModel(adata_ref)
torch.set_num_threads(os.cpu_count())

mod.train(
    max_epochs=250,
    accelerator='gpu',
    batch_size=4096,
    datasplitter_kwargs={"num_workers": 0, "persistent_workers": False},
)

# ── Export posterior ───────────────────────────────────────────────────
adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator': 'cpu'}  # cpu for export avoids VRAM OOM
)

# ── Save model + adata immediately ────────────────────────────────────
mod.save(ref_run_name, overwrite=True)

# explicitly save adata so we can recover without retraining
adata_ref.var = adata_ref.var.drop(columns=['_index'], errors='ignore')
adata_ref.obs = adata_ref.obs.drop(columns=['_index'], errors='ignore')
adata_ref.write_h5ad(f"{ref_run_name}/adata_ref_trained.h5ad")
print("Model and adata saved.")

# ── QC plots ───────────────────────────────────────────────────────────
mod.plot_history(20)
plt.savefig(f"{ref_run_name}/train_history.png", bbox_inches='tight')
plt.close('all')

mod.plot_QC()
plt.savefig(f"{ref_run_name}/train_qc.png", bbox_inches='tight')
plt.close('all')

# ── Extract signature matrix ───────────────────────────────────────────
inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][
    [f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]
].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.to_csv(f"{ref_run_name}/reference_signatures.csv")

print("Cell types found:", inf_aver.columns.tolist())
print("Signature matrix shape:", inf_aver.shape)
print(f"Done. All outputs saved to: {ref_run_name}")