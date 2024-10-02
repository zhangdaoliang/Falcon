#-*- coding : utf-8 -*-
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score
from Falcon.utils import *
from Falcon.process import *
from Falcon import train_model
from datetime import datetime

section_id = "Prostate_FFPA"
k=7

im_re = pd.read_csv("Data/Visium_FFPE_Human_Prostate_Cancer/image_representation/ViT_pca_representation.csv", header=0, index_col=0, sep=',')
print(section_id, k)
adata = sc.read_visium("Data/Visium_FFPE_Human_Prostate_Cancer",
                       count_file="Visium_FFPE_Human_Prostate_Cancer_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
adata.obsm["im_re"] = im_re

Ann_df = pd.read_csv("Data/Visium_FFPE_Human_Prostate_Cancer/Pathology.csv", sep=",", header=0,
                     index_col=0)
adata.obs['Ground Truth'] = Ann_df['Pathology']
# adata.obs['ground_truth'] = adata.obs['Ground Truth']
adata =  adata[:, adata.var['highly_variable']]
adata.obsm["adj"] = calculate_adj_matrix(adata)
adata= train_model.train(adata,k,n_epochs=25,h=[3000,3000],radius=0,l=2,
                         lr=1e-5, weight_coef=0.01, weight_selfExp=0.00001)

obs_df = adata.obs.dropna()
ARI = adjusted_rand_score(obs_df['stMGSC'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)

