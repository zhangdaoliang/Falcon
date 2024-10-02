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
from sklearn.metrics import normalized_mutual_info_score

adata = sc.read_h5ad("/home/dell/stproject/stMGSC/HD20.h5ad")
adata.var_names_make_unique()

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

adata.obsm["adj"] = calculate_adj_matrix(adata)

adata = train_model.train(adata,8,n_epochs=100,h=[1000,1000],radius=0,l=5,
                             lr=1e-6, weight_coef=0.05, weight_selfExp=0.00001)
adata.obs['stMGSC'] = adata.obs['stMGSC'].astype(str)
sc.pl.spatial(adata, color="stMGSC",save= "HD_stMGSC",title='Falcon',s=10)

# adata = sc.read_visium("Data/HD", count_file="filtered_feature_bc_matrix.h5")
# adata.var_names_make_unique()
# spatial=adata.obsm["spatial"]
# used_barcode=pd.read_csv("Data/HD/cell_used.csv",header=None,index_col=0)
# used_barcode = used_barcode.index
# adata = adata[used_barcode,]
# im_re = pd.read_csv('Data/HD/images_representation/ViT_pca_representation.csv',
#                     header=0, index_col=0, sep=',')
# adata.obsm["im_re"] = im_re
#
# # 设定你要删除的细胞比例
# fraction_to_remove = 0.8  # 例如删除10%的细胞
#
# # 确定需要删除的细胞数量
# num_cells_to_remove = int(adata.n_obs * fraction_to_remove)
#
# # 随机选择要删除的细胞索引
# # np.random.seed(42)  # 设定随机种子以获得可重复的结果
# cells_to_remove = np.random.choice(adata.obs_names, num_cells_to_remove, replace=False)
#
# # 删除这些细胞
# adataT = adata[~adata.obs_names.isin(cells_to_remove)].copy()
#
# sc.pp.highly_variable_genes(adataT, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_per_cell(adataT)
# sc.pp.log1p(adataT)
#
# adataT.obsm["adj"] = calculate_adj_matrix(adataT)
#
# adataT = train_model.train(adataT,10,n_epochs=3,h=[3000,3000],radius=0,l=0.5,
#                              lr=1e-6, weight_coef=0.05, weight_selfExp=0.00001)
# adataT.obs['stMGSC'] = adataT.obs['stMGSC'].astype(str)
# sc.pl.embedding(adataT , basis="spatial", color="stMGSC", s=6, show=True,save = "HD_stMGSC")

