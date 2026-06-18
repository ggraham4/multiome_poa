import scanpy as sc
import scipy.io
import pandas as pd

counts = scipy.io.mmread("C:/Users/Gabe/Desktop/Visium/counts.mtx").T.tocsr()
obs    = pd.read_csv("C:/Users/Gabe/Desktop/Visium/meta.csv", index_col=0)
genes  = pd.read_csv("C:/Users/Gabe/Desktop/Visium/genes.csv", index_col=0)

# genes.csv has a column called 'gene' — use that as the index
var = pd.DataFrame(index=genes["gene"].values)

adata_vis = sc.AnnData(X=counts, obs=obs, var=var)

print(adata_vis)
print("var_names sample:", adata_vis.var_names[:5].tolist())

# Load reference to check overlap
import pandas as pd
inf_aver = pd.read_csv("C:/Users/Gabe/Desktop/multiome_poa/Python/cell2location/results_2026_06_17/reference_signatures/reference_signatures.csv", index_col=0)
shared = adata_vis.var_names.intersection(inf_aver.index)
print(f"Shared genes: {len(shared)}")

adata_vis.write_h5ad("C:/Users/Gabe/Desktop/Visium/vis_dissection.h5ad")
print("Done.")