# CellTypist benchmark: GSE155960 adipose macrophage-enriched subset
# CellTypist 1.7.1; Immune_All_Low.pkl

from pathlib import Path

import celltypist
import pandas as pd
import scanpy as sc
import scipy.io
import numpy as np
from celltypist import models

base_dir = Path("benchmark_CellTypist/GSE155960")
input_dir = base_dir / "CellTypist_input"
input_dir.mkdir(parents=True, exist_ok=True)

# Load the raw gene-by-cell count matrix exported from Seurat.
matrix = scipy.io.mmread(input_dir / "matrix.mtx").tocsr()
genes = pd.read_csv(input_dir / "genes.tsv", header=None)[0].astype(str)
cells = pd.read_csv(input_dir / "cells.tsv", header=None)[0].astype(str)

adata = sc.AnnData(matrix.T)
adata.var_names = genes
adata.obs_names = cells
adata.layers["counts"] = adata.X.copy()

# Apply CellTypist preprocessing.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Annotate cells using the Immune_All_Low reference model.
model = models.Model.load(model="Immune_All_Low.pkl")
predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=False
)

# Save predicted labels and maximum class score.
labels = predictions.predicted_labels.iloc[:, 0]
class_scores = predictions.probability_matrix

results = pd.DataFrame({
    "Cell": labels.index,
    "CellTypist_label": labels.values,
    "CellTypist_probability": class_scores.max(axis=1).reindex(labels.index).values
})
results.to_csv(base_dir / "CellTypist_predictions.csv", index=False)

# Save the complete cell-by-class score matrix.
class_scores.to_csv(base_dir / "CellTypist_probabilities.csv", index_label="Cell")
