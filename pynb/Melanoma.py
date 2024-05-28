# %%
# Import necessary libraries
import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix, lil_matrix
from pathlib import Path

inputfile = "../data/Melanoma/GSE72056_melanoma_single_cell_revised_v2.txt"
outputfile = "../out/Melanoma"

HVF = True
TEST = True

TEST_genes = 800
TEST_samples = 10000

k = 15

# %% [markdown]
# # Load Dataset

# %%
se = pd.read_csv(inputfile, index_col=0, sep="\t")[3:]
meta = pd.read_csv(inputfile, index_col=0, sep="\t").iloc[:3]
se = csr_matrix(se.values)
se = se[se.getnnz(1) > 0][:, se.getnnz(0) > 0]

if TEST:
    tgenes = min(TEST_genes, se.shape[1])
    tsamples = min(TEST_samples, se.shape[0])
    se.resize((tsamples, tgenes))
    se = se[se.getnnz(1) > 0][:, se.getnnz(0) > 0]

adata = sc.AnnData(se)
# sc.pp.normalize_total(adata, target_sum=1e6)
# sc.pp.log1p(adata)
sc.pp.scale(adata, zero_center=False)
sc.pp.highly_variable_genes(adata)


# %%
if HVF:
    adata = adata[:, adata.var.highly_variable]
    # adata = adata[csr_matrix(adata.X).getnnz(1) > 0][:, csr_matrix(adata.X).getnnz(0) > 0]
    adata = adata[adata.X.getnnz(1) > 0][:, adata.X.getnnz(0) > 0]

print(adata)

sc.pp.neighbors(adata, n_neighbors=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# %% [markdown]
# # Visualize data

# %%
# PCA plot
imgname = os.path.join(outputfile, "PCA.png")
print(f"Saving Image --- {imgname}")
sc.pl.pca(adata, save=imgname.split("/")[-1])

# %%
# UMAP plot
imgname = os.path.join(outputfile, "UMAP.png")
print(f"Saving Image --- {imgname}")
sc.pl.umap(adata, save=imgname.split("/")[-1])

# %%
# Elbow plot
imgname = os.path.join(outputfile, "elbow.pdf")
print(f"Saving Image --- {imgname}")
sc.pl.pca_variance_ratio(adata, log=True, save=imgname.split("/")[-1])

# %% [markdown]
# # Archetypes

# %%
import archetypes as arch
from time import time

aa_kwargs = {
    "n_archetypes": 4,
    "n_init": 5,
    "max_iter": 100000,
    "verbose": True,
    "tol": 1e-3,
}

mod0 = arch.AA(**aa_kwargs, algorithm_init="furthest_sum")

t0 = time()
mod0.fit_transform(adata.X.toarray())
t1 = time()

print(f"mod0: {t1-t0:.2f} seconds|RSS: {mod0.rss:.2f}")


# %%
