# runtime benchmarking wrt k-means
# nx-metis installation is recommended via colab
# for which the notebook is provideds

import os
import sys
from pathlib import Path
import time
import random
import torch
import numpy as np
import scipy.io as sp
import networkx as nx
from sklearn.utils import shuffle
from sklearn.cluster import KMeans
from torch_geometric.data import Data
import matplotlib.pyplot as plt

sys.path.append('../')       

from models import SageBase, SageRes
from losses import loss_normalized_cut
from utils import get_dataset, get_sample

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

adjacencies, coords, areas, dataset_size = get_dataset('timebenchmark',True)
size = np.array([len(coords[i]) for i in range(len(coords))])
ntests = 20

sagebase_time = torch.zeros((ntests,dataset_size))

model = SageBase(64,32,3,2)
model.load_state_dict(torch.load('../outputs/models/SageBase.pt'))

for j in range(ntests):
  for i in range(dataset_size):
    data_test = get_sample(adjacencies, coords, areas,i)
    start = time.monotonic()
    out = model(data_test.x,data_test.edge_index)
    sagebase_time[j,i] = time.monotonic()-start


sageres_time = torch.zeros((ntests,dataset_size))

model = SageRes()
model.load_state_dict(torch.load('../outputs/models/SageRes.pt'))
model.eval()

for j in range(ntests):
  for i in range(dataset_size):
    data_test = get_sample(adjacencies, coords, areas,i)
    start = time.monotonic()
    out = model(data_test.x,data_test.edge_index)
    sageres_time[j,i] = time.monotonic()-start


kmeans_time = torch.zeros((ntests,dataset_size))

for j in range(ntests):
  for i in range(dataset_size):
    start = time.monotonic()
    out = KMeans(n_clusters=2, random_state=0).fit(coords[i])
    kmeans_time[j,i]= time.monotonic()-start


sagebase_time_avg = torch.mean(sagebase_time,0)
sageres_time_avg = torch.mean(sageres_time,0)
#metis_time_avg = torch.mean(metis_time,0)
kmeans_time_avg = torch.mean(kmeans_time,0)

m1, b1 = np.polyfit(size, sagebase_time_avg, 1)
m2, b2 = np.polyfit(size, sageres_time_avg, 1)
#m3, b3 = np.polyfit(size, metis_time_avg, 1)
m4, b4 = np.polyfit(size, kmeans_time_avg, 1)

plt.figure(figsize=(10,6))
plt.loglog(size,size*m1+b1,'-',color='#2CDA9D')
plt.loglog(size,size*m2+b2,'-',color='#F26DF9')
#plt.semilogy(size,size*m3+b3,'-',color='#2F97C1')
plt.loglog(size,size*m4+b4,'-',color='#FF1053')

for i in range(ntests):
  plt.loglog(size,sagebase_time[i,:],'o',color='#2CDA9D',alpha=.2)
  plt.loglog(size,sageres_time[i,:],'o',color='#F26DF9',alpha=.2)
 # plt.semilogy(size,metis_time[i,:],'o',color='#2F97C1',alpha=.2)
  plt.loglog(size,kmeans_time[i,:],'o',color='#FF1053',alpha=.2)


fs = 16

plt.legend(['SAGE-Base','SAGE-Res','k-means'],fontsize=fs) #'METIS',
plt.xlabel('Graph Size',fontsize=fs)
plt.ylabel('Time [s]',fontsize=fs)
plt.title('Bisection Methods Runtime',fontsize=fs)
plt.savefig('runtime.png',dpi=300)