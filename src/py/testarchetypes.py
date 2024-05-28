# %%
import pandas as pd
import scanpy as sc


def loadMelanoma(hvf_nfeatures=2000):
    filename = "/app/data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt"
    metadata = pd.read_table(filename, header=0, index_col=0, nrows=3, sep="\t")
    metadata = metadata.T.apply(pd.to_numeric, errors="coerce")

    # Basic setup
    ge = pd.read_table(filename, header=0, index_col=0, skiprows=[1, 2, 3], sep="\t")
    print(f"size pre deduplicate {ge.shape}")

    ge = ge[~ge.index.duplicated(keep="first")]
    print(f"size after deduplicate {ge.shape}")

    # Seurat
    adata = sc.AnnData(X=ge.values.T, obs=metadata)
    # sc.pp.scale(adata, zero_center=True, max_value=10)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvf_nfeatures)
    sc.pp.pca(adata)

    return adata


adata = loadMelanoma()
adata_hvf = adata[:, adata.var["highly_variable"].values]

# %%
import archetypes as arch

aa_kwargs = {
    "n_archetypes": 20,
    "n_init": 1,
    "max_iter": 100000,
    "verbose": True,
    "tol": 1e-4,
}
# model= arch.AA(**aa_kwargs, algorithm_init="random")
model = arch.AA(**aa_kwargs, algorithm_init="furthest_sum")

# %%
model.fit_transform(adata_hvf.X)

model.archetypes_
