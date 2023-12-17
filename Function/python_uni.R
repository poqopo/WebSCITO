# 필요한 라이브러리 불러오기

library(data.table)
library(reticulate)
use_python("/home/poqopo/anaconda3/bin/python")

py_run_string("
import re
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
import leidenalg
sc.settings.verbosity = 4

adt_2_df = r.selected_ab_cell_matrix
meta_2_df = r.meta_ab_cell_matrix
adt_2_df = adt_2_df.astype('float')
adt_adata = sc.AnnData(adt_2_df.T)
adt_adata.obs = meta_2_df

sc.pp.normalize_per_cell(adt_adata, counts_per_cell_after=1e4)

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split('(\\d+)', text)]


def normalize_byBatch(adt, batch_var='Batch', max_val=10):
    all_batches = list(set(adt.obs[batch_var]))
    all_batches.sort(key=natural_keys)
    ann_dats = []
    for b in all_batches:
        batch_adat = adt[adt.obs[batch_var] == b]
        sc.pp.log1p(batch_adat)
        sc.pp.scale(batch_adat, max_value=max_val)
        ann_dats.append(batch_adat)
    norm_concat = ann_dats[0].concatenate(ann_dats[1:len(ann_dats)+1])
    return norm_concat

adt_adata = normalize_byBatch(adt_adata)

adt_adata.var_names_make_unique()


ptn_names =r.ptn_names

adt_adata.var = ptn_names


sc.tl.pca(adt_adata, svd_solver='arpack', use_highly_variable=None)
sc.pp.neighbors(adt_adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adt_adata)
sc.tl.leiden(adt_adata)
")
