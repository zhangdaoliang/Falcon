# -*- coding: utf-8 -*-
import scipy.sparse as sp
import torch
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.cluster import KMeans
from Falcon.process import *
from Falcon import models
from Falcon.loss_function import loss_function
import numpy as np
import torch.optim as optim
from Falcon.utils import *
import warnings
warnings.filterwarnings('ignore')
from sklearn.decomposition import PCA
from tqdm import tqdm
def train(adata,knn=10,h=[1000,1000], n_epochs=200,lr=0.000001, key_added='stMGSC', random_seed=0,
          l=2,weight_decay=0.0001, weight_coef=0.01, weight_selfExp=0.00001,radius=0,
                device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')):
    set_seed(random_seed)
    if 'highly_variable' in adata.var.columns:
        adata_Vars =  adata[:, adata.var['highly_variable']]
    # labels = adata_Vars.obs['ground_truth']
    features_X = torch.FloatTensor(adata_Vars.X.toarray()).to(device)
    features_I = torch.FloatTensor(adata_Vars.obsm["im_re"].values).to(device)
    adj=adata_Vars.obsm["adj"]
    adj = np.exp(-1*(adj**2)/(2*(l**2)))
    adj = sp.coo_matrix(adj)
    adj = normalize(adj + sp.eye(adj.shape[0]))
    adj = sparse_mx_to_torch_sparse_tensor(adj).to(device)
    model =models.stMGSC(nfeatX=features_X.shape[1],
                 nfeatI=features_I.shape[1],
                    hidden_dims=h,
                    num_sample=features_X.shape[0]
                    ).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    # optimizer = optim.SGD(model.parameters(), lr=lr, weight_decay=weight_decay)

    for epoch in tqdm(range(n_epochs), desc="Training Progress"):
        model.train()
        optimizer.zero_grad()
        f_list, zf_list = model(features_X, features_I, adj)
        loss = loss_function(model.self_expression.Coefficient, f_list, zf_list, weight_coef, weight_selfExp)
        loss.backward()
        optimizer.step()

    z_xi = torch.cat(f_list, dim=1)

    kmeans = KMeans(n_clusters=knn,random_state=random_seed).fit(np.nan_to_num(z_xi.cpu().detach()))
    idx_max = kmeans.labels_
    emb_max = z_xi.to('cpu').detach().numpy()

    adata.obs["cluster"] = idx_max.astype(str)
    if radius !=0 :
        nearest_new_type = refine_label(adata, radius=radius)
        adata.obs[key_added] = nearest_new_type
    else:
        adata.obs[key_added] = adata.obs["cluster"]
    # adata.obsm["emb"] = emb_max
    pca = PCA(n_components=50, random_state=random_seed)
    adata.obsm['emb_pca'] = pca.fit_transform(emb_max.copy())
    return adata

