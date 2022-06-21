import pyVIA.core as via
from scipy.io import loadmat
import pandas as pd
import scanpy as sc
import numpy as np
import warnings

adata=sc.read_h5ad()

time_labels=

input_data = adata.obsm['X_pca'][:, 0:ncomps]
v0 = via.VIA(input_data, time_labels, jac_std_global=0.15, dist_std_local=1, knn=knn,
             cluster_graph_pruning_std=.15,
             too_big_factor=v0_too_big, root_user=root_user, dataset='EB', random_seed=v0_random_seed,
             do_impute_bool=True, is_coarse=True, preserve_disconnected=True) 
v0.run_VIA()
