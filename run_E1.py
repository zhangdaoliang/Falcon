#-*- coding : utf-8 -*-
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score
from Falcon.utils import *
from Falcon.process import *
from Falcon import train_model
from datetime import datetime
# for d in ["A1","B1","C1","D1","E1","F1","G2","H1"]:
section_id = "A1"
# k=5
# l=100
# 0.45
section_id = "B1"
# k=4
# l=25000
# # 0.28
section_id = "C1"
# k=3
# l=25000
section_id = "D1"
# k=3
# l=20000
section_id = "E1"
# k=3
# l=200
section_id = "F1"
# k=3
# l=500
# section_id = "G2"
k=6
l=200
section_id = "H1"
# k=6
# l=200
# for section_id, k ,l in zip(["A1","B1","C1","D1","E1","F1","G2","H1"], [5,4,3,3,3,3,6,6],
#                             [100,25000,25000,20000,200,500,200,200]):

im_re = pd.read_csv("Data/STbreast_Tumors/{}/image_representation/ViT_pca_representation.csv".format(section_id),
                    header=0, index_col=0, sep=',')
print(section_id, k)
genedata= pd.read_csv("Data/STbreast_Tumors/{}/ut_{}_stdata_filtered.tsv".format(section_id,section_id),
                      header=0, sep='	',index_col=0)
labeldf = pd.read_csv("Data/STbreast_Tumors/{}/{}_label.csv".format(section_id,section_id),
                      sep=",", header=0,na_filter=False, index_col=0)
labeldf.replace("undetermined", pd.NA, inplace=True)
adata = sc.AnnData(genedata)
adata.var_names_make_unique()
prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
prefilter_specialgenes(adata)
#Normalization
# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# sc.pp.scale(adata, zero_center=False, max_value=10)

adata.obs['Ground Truth'] = labeldf["label"]
adata.obs['ground_truth'] = adata.obs['Ground Truth']
prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_per_cell(adata)
# sc.pp.log1p(adata)
adata.obsm["im_re"] = im_re
adata.obs["y"] = labeldf["pixel_y_y"]
adata.obs["x"] = labeldf["pixel_x_y"]
adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
adata =  adata[:, adata.var['highly_variable']]
a=adata.obsm["spatial"]
adata.obsm["adj"] = calculate_adj_matrix(adata)
adata= train_model.train(adata,k,n_epochs=500,h=[3000,3000],radius=0,l=l,
                         lr=1e-6, weight_coef=0.5, weight_selfExp=0.0005)
obs_df = adata.obs.dropna()
# obs_df.to_csv("result/{}_type_stMMR.csv".format(section_id))
ARI = adjusted_rand_score(obs_df['stMGSC'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)

