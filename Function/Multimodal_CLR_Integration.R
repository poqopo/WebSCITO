library(Matrix)
library(Seurat)
library(cluster)
library(fitdistrplus)
library(data.table)
library(dplyr)
library(ggrepel)
library(dsb)
options(future.globals.maxSize = 8000 * 1024^2)
library(reticulate)
use_python("/home/poqopo/anaconda3/bin/python")

raw <- readRDS("/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/raw.rds")
ptn_names <- read.table("/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/ptn_names.csv")
meta <- read.table("/data/project/scitoseq/script/Final/Sample/Multi_4AB_4Pool/batch_meta.csv")
abNum <- 209
batchNum <- 4

  tx <- raw$`Gene Expression`[, which(Matrix::colSums(raw$`Gene Expression`) > 0)]
  ab <- raw$`Antibody Capture`[, which(Matrix::colSums(raw$`Antibody Capture`) > 0)]
  ccd_indices <- intersect(which(Matrix::colSums(tx) > 0), which(Matrix::colSums(ab) > 0))
  
  ab <- ab[, ccd_indices]
  tx <- tx[, ccd_indices]
  
  ab_resolved_list <- list()
  tx_resolved_list <- list()
  
  for (i in 1:as.integer(batchNum)) {
    start <- (1 + abNum*(i-1))
    end <- i*abNum
    
    ab_subset <- ab[start:end, which(meta$Batch == paste("Batch",i,sep=""))]
    ab_resolved_list[[i]] <- ab_subset
    
    tx_subset <- tx[, which(meta$Batch == paste("Batch",i,sep=""))]
    tx_resolved_list[[i]] <- tx_subset
  }
  
  ab_resolved <- do.call(cbind, ab_resolved_list)
  tx_resolved <- do.call(cbind, tx_resolved_list)

  cell_names <-make.names(colnames(ab_resolved), unique = TRUE)
  gene_names <- make.names(rownames(tx_resolved), unique = TRUE)
  ptn_names <- ptn_names$Antibody
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

tx_df = r.tx_resolved
tx_df = tx_df.astype('float')
rna = sc.AnnData(tx_df.T)
ab_df = r.ab_resolved
ab_df = ab_df.astype('float')
prot = sc.AnnData(ab_df.T)
rna.obs_names= r.cell_names
rna.var_names = r.gene_names
prot.obs_names = r.cell_names
prot.var_names = r.ptn_names
rna.obs['Pool'] = r.meta

og_shape = rna.shape

sc.pp.filter_cells(rna, min_counts=1, inplace=True)
rem_cell_shape = rna.shape

sc.pp.filter_genes(rna, min_counts=1, inplace=True)
rem_cell_gene_shape = rna.shape
rna.var['mt'] = rna.var_names.str.startswith('MT-')
rna.var['rb'] = rna.var_names.str.startswith(('RPS', 'RPL'))
rna.var['hb'] = rna.var_names.str.contains(('^HB[^(P)]'))
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=['mt', 'rb', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)
rna.obs['complexity'] = np.log10(rna.obs['n_genes_by_counts']) / np.log10(rna.obs['total_counts'])
rna.obs['complexity'] = rna.obs['complexity'].astype(np.float32)

fil_dict = dict(
    total_counts=(500, 30000),
    n_genes_by_counts=(500, 6000),
    complexity=(0.88,),
    pct_counts_mt=(7.5,),
    pct_counts_hb=(0.2,),
    gene_thres=20
)
print('Performing cell-level filtering...')

mu.pp.filter_obs(rna, 'total_counts', lambda x: x >= fil_dict['total_counts'][0])
mu.pp.filter_obs(rna, 'total_counts', lambda x: x < fil_dict['total_counts'][1])
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: x >= fil_dict['n_genes_by_counts'][0])
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: x < fil_dict['n_genes_by_counts'][1])
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < fil_dict['pct_counts_mt'])
mu.pp.filter_obs(rna, 'pct_counts_hb', lambda x: x < fil_dict['pct_counts_hb'])
gene_cell_thres = int(fil_dict['gene_thres'])
sc.pp.filter_genes(rna, min_cells=gene_cell_thres)
sc.pp.calculate_qc_metrics(
    rna,
    qc_vars=['mt', 'rb', 'hb'],
    percent_top=None,
    log1p=False,
    inplace=True
)
rna.obs['complexity'] = np.log10(rna.obs['n_genes_by_counts']) / np.log10(rna.obs['total_counts'])
cols = ['total_counts', 'total_counts_mt', 'total_counts_rb', 'total_counts_hb']
rna.obs[cols] = rna.obs[cols].astype(np.int32)
rna.obs['complexity'] = rna.obs['complexity'].astype(np.float32)
cols = ['total_counts', 'n_cells_by_counts']
rna.var[cols] = rna.var[cols].astype(np.int32)
cols = ['pct_dropout_by_counts']
rna.var[cols] = rna.var[cols].astype(np.float32)
for x in ('n_cells', 'n_counts'):
    if x in rna.var_keys():
        del(rna.var[x])
    if x in rna.obs_keys():
        del(rna.obs[x])
mdata = MuData({'rna': rna, 'prot': prot})
mdata.mod['rna'] = rna
mu.pp.intersect_obs(mdata)
  ")
py_run_string("
sc.pp.filter_genes(prot, min_counts=1, inplace=True)
rem_cell_gene_shape = prot.shape
pt.pp.clr(prot)
prot.layers['clr_counts'] = prot.X.copy()
prot.obs['total_clr_counts'] = prot.X.sum(axis=1)
mdata.mod['prot'] = prot
rna.layers['counts'] = rna.X.copy()
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
rna.layers['logcounts'] = rna.X.copy()

print('Computing PCA...')
sc.pp.pca(rna)

rna.layers['scaled'] = rna.X.copy()
rna.X = rna.layers['logcounts'].copy()
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=15)
sc.tl.umap(rna)
sc.tl.leiden(rna)
mdata.mod['rna'] = rna
rna = mdata['rna'].copy()
prot = mdata['prot'].copy()
")

py_run_string("
mu.pl.embedding(
    mdata,
    basis='rna:umap',
    color=['CD4','CD4_protein'],
    cmap=cc.cm.gouldian,
    palette=cc.glasbey_hv)
")
