import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
from sklearn.cluster import KMeans
import torch.optim as optim
from random import shuffle
import pandas as pd
import numpy as np
import scanpy as sc
from Falcon.layers import GraphConvolution,SelfAttention,MLP,MGCN,decoder
import sys
class SelfExpression(nn.Module):
    def __init__(self, n, weight_c):
        super(SelfExpression, self).__init__()
        self.Coefficient = nn.Parameter(weight_c * torch.randn(n, n), requires_grad=True)

    def forward(self, x):  # shape=[n, d]
        y = torch.matmul(self.Coefficient, x)
        return y

class stMGSC(nn.Module):
    def __init__(self,nfeatX,nfeatI,hidden_dims,num_sample, weight_c = 1e-2):
        super(stMGSC, self).__init__()
        self.n = num_sample
        self.mgcn = MGCN(nfeatX,nfeatI,hidden_dims)
        self.self_expression = SelfExpression(self.n, weight_c)

    def forward(self,x,i,a):
        emb_x,emb_i = self.mgcn(x,i,a)
        zx = self.self_expression(emb_x)
        zi = self.self_expression(emb_i)
        return [emb_x,emb_i],[zx,zi]


