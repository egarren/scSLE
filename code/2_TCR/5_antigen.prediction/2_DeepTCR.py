import matplotlib
matplotlib.use("agg") 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from DeepTCR.DeepTCR import DeepTCR_U
from DeepTCR.DeepTCR import DeepTCR_SS

###Unsupervised CNN on repertoires (disease-based) for clustering
# Instantiate training object
DTCRUr = DeepTCR_U('./unsupervised.rep')
#Load Data from directories
DTCRUr.Get_Data(directory='../SLE.data',Load_Prev_Data=False,aggregate_by_aa=True,
               aa_column_beta=2,v_beta_column=3,j_beta_column=4, d_beta_column=5,
               aa_column_alpha=6, v_alpha_column=7,j_alpha_column=8,count_column=12)
colors={"B6":"grey","SLE.yaa":"red"}
# #Train VAE
DTCRUr.Train_VAE(Load_Prev_Data=False,var_explained=0.99)
features = DTCRUr.features
DTCRUr.Sample_Features()
DTCRUr.sample_features
print(features.shape)
DTCRUr.sample_features.to_csv("unsupervised.rep.features.csv")
DTCRUr.HeatMap_Sequences()
plt.savefig("unsupervised.rep.heatmap.seq.png")
DTCRUr.HeatMap_Samples()
plt.savefig("unsupervised.rep.heatmap.sample.ratio.png")
#Clustering
DTCRUr.Cluster(clustering_method='phenograph',write_to_sheets=True,sample=1000) #hierarchical, dbscan; use sample to downsample
DFs = DTCRUr.Cluster_DFs
df = pd.concat(DFs,keys=list(range(0,len(DFs))),names=['cluster','rowID'])
df.to_csv("unsupervised.rep.clusters.csv")
DTCRUr.UMAP_Plot(by_cluster=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.clusters.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_sample=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.samples.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5)
plt.savefig("unsupervised.rep.umap.VAE.class.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5,show_legend=False)
plt.savefig("unsupervised.rep.umap2.VAE.class.png")
DTCRUr.UMAP_Plot_Samples(scale=100)
plt.savefig("unsupervised.rep.umap.samples.png")
#export UMAP
import umap
import numpy
import pickle
umap_obj = umap.UMAP()
features_umap = umap_obj.fit_transform(DTCRUr.features)
df=numpy.c_[features_umap,DTCRUr.class_id,DTCRUr.sample_id,DTCRUr.beta_sequences,DTCRUr.alpha_sequences,DTCRUr.freq,DTCRUr.counts]
numpy.savetxt("unsup.rep.umap.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
pickle.dump(DTCRUr.features,open("unsup.rep.obj","wb"))
#Repertoire visualization
DTCRUr.Repertoire_Dendrogram(n_jobs=40,distance_metric='KL')
DTCRUr.Repertoire_Dendrogram(lw=6,gridsize=25,Load_Prev_Data=True,
    gaussian_sigma=0.75,dendrogram_radius=0.2,repertoire_radius=0.3,color_dict=colors) #,sample_labels=True
plt.savefig("unsupervised.rep.repertoire.dendrogram.png")

