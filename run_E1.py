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
# 78 epoch=78
obs_df = adata.obs.dropna()
# obs_df.to_csv("result/{}_type_stMMR.csv".format(section_id))
ARI = adjusted_rand_score(obs_df['stMGSC'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)
# now = datetime.now()
# dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
# with open("Data/out.txt", "a") as file:
#     file.write(
#         dt_string + "---,dataset:{},ARI:{},r:{}\n"
#         .format(section_id, ARI,r)
#         )

#
ax=sc.pl.scatter(adata, alpha=1, x="x", y="y", color="stMGSC", legend_fontsize=18, show=False,
                   size=50000 / adata.shape[0])

title='stMGSC (ARI=%.2f)' % ARI
ax.set_title(title, fontsize=23)
ax.set_aspect('equal', 'box')
ax.set_xticks([])
ax.set_yticks([])
ax.axes.invert_yaxis()
plt.xlabel('')
plt.ylabel('')
plt.savefig("{}_{}.pdf".format(section_id, "stMGSC"))
plt.show()


ax=sc.pl.scatter(adata, alpha=1, x="x", y="y", color="Ground Truth", legend_fontsize=18, show=False,title='Ground Truth',
                   size=50000 / adata.shape[0])
title='Ground Truth'
ax.set_title(title, fontsize=23)
ax.set_aspect('equal', 'box')
ax.set_xticks([])
ax.set_yticks([])
ax.axes.invert_yaxis()
plt.xlabel('')
plt.ylabel('')
plt.savefig("{}_{}.pdf".format(section_id, "GT"))
plt.show()




# sc.pl.scatter(adata, alpha=1, x="x", y="y", color=["stMGSC", "Ground Truth"], legend_fontsize=18, show=True,
#     title=['stMGSC (ARI=%.2f)' % ARI, "Ground Truth"],
#                    size=100000 / adata.shape[0],save=section_id)

# plt.rcParams["figure.figsize"] = (3, 3)
# sc.pl.spatial(adata, color=["stMGSC", "Ground Truth"], title=['stMGSC (ARI=%.2f)' % ARI, "Ground Truth"],
#            save=section_id)
# #
# sc.pp.neighbors(adata, use_rep='emb_pca')
# sc.tl.umap(adata)
# plt.rcParams["figure.figsize"] = (3, 3)
# sc.pl.umap(adata, color=["stMGSC", "Ground Truth"], title=['stMGSC (ARI=%.2f)' % ARI, "Ground Truth"],
#            save=section_id,size=100)
#
# sc.tl.rank_genes_groups(adata, "stMGSC", method='wilcoxon',corr_method="benjamini-hochberg")

# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# df  = sc.get.rank_genes_groups_df(adata, group="4",pval_cutoff=0.01,log2fc_min=1.5)
# df = df.sort_values(by="scores", ascending=False)
# df.to_csv('HBD4.csv')
# print(df)
