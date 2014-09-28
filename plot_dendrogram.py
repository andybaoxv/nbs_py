import numpy as np
from python.COPDGene.utils.my_hierarchical_clustering import\
        my_hierarchical_clustering
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import copy

basepath = "/home/changyale/matlab/nbs_release_v0.2/result/"
filename = "NMF_cooccurence_matrix.csv"
#file_label = "NMF_ID_labels.csv"

mtr_sim = np.loadtxt(basepath+filename,delimiter=',')
mtr_lin = my_hierarchical_clustering(mtr_sim,method='average')

# Draw dendrogram of hierarchical clustering
fig = plt.figure(figsize=(30,30))
tmp = copy.deepcopy(mtr_lin)
for i in range(tmp.shape[0]):
    tmp[i,2] = np.max(mtr_dis)-tmp[i,2]
dend = dendrogram(tmp,orientation='right',line_width=10)
plt.savefig('dendrogram_hierarchical_ecoli.png')

