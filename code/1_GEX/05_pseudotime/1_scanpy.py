import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import rpy2.robjects as robjects
import anndata2ri

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/scanpy.h5ad'
anndata2ri.activate()
robjects.r['load']('sce.RData')
adata = robjects.r('as(SLE.B.obj.combined.sce, "SingleCellExperiment")')
adata = robjects.r('as(SLE.CD4.obj.combined.sce, "SingleCellExperiment")')
adata = robjects.r('as(SLE.CD8.obj.combined.sce, "SingleCellExperiment")')
#compute PAGA path 
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
sc.tl.paga(adata, groups='cluster_letter')
sc.pl.paga(adata,save=True)
#plot PAGA path
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    palette=["#F8766D", "#CD9600" ,"#7CAE00", "#00BE67", "#00BFC4" ,"#00A9FF","#C77CFF","#FF61CC"], #",
    legend_fontsize=0, fontsize=0, frameon=False, edges=True, save=True)

#Diffmap
sc.tl.diffmap(adata)
diffmap = adata.obsm['X_diffmap']
np.savetxt("scanpy.diffmap.txt", diffmap, '%5.10f',delimiter = "\t")

#plot pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['cluster_letter']  == 'A')[0]
sc.tl.dpt(adata)
dpt = adata.obs.dpt_pseudotime
dpt.to_csv(path="scanpy.dpt.txt", sep="\t", float_format = '%5.10f')

