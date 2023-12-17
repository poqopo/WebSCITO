library(reticulate)
use_python("/home/poqopo/anaconda3/bin/python")

py_run_string("
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy import stats
import scanpy.external as sce
import muon as mu
from muon import MuData
from muon import prot as pt
import helper as _hp
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import anndata as ad
import colorcet as cc
")
py_run_string("
mdata = mu.read_h5mu('/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/umap.h5mu')
rna = mdata['rna'].copy()
prot = mdata['prot'].copy()
")

py_run_string("print(prot.var_names)")
py_run_string("
prot.obs['leiden'] = rna.obs['leiden']
ncols=2
nrows=1
figsize=4
wspace=0.5
fig,axs = plt.subplots(nrows=nrows, ncols=ncols,
                       figsize=(ncols*figsize+figsize*wspace*(ncols-1),nrows*figsize))
plt.subplots_adjust(wspace=wspace)
sc.pl.dotplot(rna, groupby = 'leiden', var_names = ['CD19', 'CD4'], show =False, ax=axs[0])
sc.pl.dotplot(prot, groupby = 'leiden', var_names = ['CD45_protein', 'CD4_protein'], ax=axs[1])
")
py_run_string(paste0("sc.pl.umap(adt_adata, color=['",input$selected_umap_option,"'], show = False)"))
py_run_string("plt.subplots_adjust(right=0.7)")
py_run_string("
mu.pl.embedding(
    mdata,
    basis='rna:umap',
    color=['CD19','CD19_protein'],
    cmap=cc.cm.gouldian,
    palette=cc.glasbey_hv)
")
