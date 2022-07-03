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
from torch_geometric.data import Data
import matplotlib.pyplot as plt

sys.path.append('../')       

from models import SageBase, SageRes
from losses import loss_normalized_cut
from utils import get_dataset, get_sample

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

batch_sizes = [4,8,16,32,64]
epoch_loss_container = []

from numpy import genfromtxt
epoch_loss_container = genfromtxt('data/loss_batchsize.csv', delimiter=',')

plt.figure(figsize=(10,5))
#for i in range(len(batch_sizes)):

plt.plot(np.array(epoch_loss_container[0])/epoch_loss_container[0][0])  
plt.savefig("plots/loss_batchsize_test_1.png",dpi=300)
plt.figure(figsize=(10,5))

plt.plot(np.array(epoch_loss_container[1])/epoch_loss_container[1][0])  

plt.title('Hyperparameters search - batch size')
plt.xlabel('Epochs')
#plt.legend(batch_sizes) 
plt.savefig("plots/loss_batchsize_test_2.png",dpi=300)