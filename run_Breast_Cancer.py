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

section_id = "V1_Breast_Cancer_Block_A_Section_1"
k=20

im_re = pd.read_csv(os.path.join('Data',section_id,
        "image_representation/VIT_pca_representation.csv"), header=0, index_col=0, sep=',')
print(section_id, k)
adata = sc.read_visium("Data/V1_Breast_Cancer_Block_A_Section_1",
                count_file="V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
adata.obsm["im_re"] = im_re

Ann_df = pd.read_csv("Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                     index_col=0)
adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']
# adata.obs['ground_truth'] = adata.obs['Ground Truth']
adata =  adata[:, adata.var['highly_variable']]
adata.obsm["adj"] = calculate_adj_matrix(adata)


adata= train_model.train(adata,k,n_epochs=20,h=[3000,3000],radius=50,l=0.76,
                         lr=1e-6, weight_coef=0.01, weight_selfExp=0.00001)
# adata.write("result/breast/breast.h5ad")
obs_df = adata.obs.dropna()
# obs_df.to_csv("result/{}_type_stMMR.csv".format(section_id))
ARI = adjusted_rand_score(obs_df['stMGSC'], obs_df['Ground Truth'])



print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['stMGSC'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, color=["stMGSC", "Ground Truth"], title=['stMGSC (ARI=%.2f)' % ARI, "Ground Truth"])
#
# sc.pp.neighbors(adata, use_rep='emb')
# sc.tl.umap(adata)
# plt.rcParams["figure.figsize"] = (3, 3)
# sc.pl.umap(adata, color=["stMGSC", "Ground Truth"], title=['stMGSC (ARI=%.2f)' % ARI, "Ground Truth"],
#            save="breast.pdf")
#
# sc.tl.rank_genes_groups(adata, "stMGSC", method='wilcoxon',corr_method="benjamini-hochberg")
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# df  = sc.get.rank_genes_groups_df(adata, group="4",pval_cutoff=0.01,log2fc_min=1.5)
# df = df.sort_values(by="scores", ascending=False)
# df.to_csv('HBD4.csv')
# print(df)