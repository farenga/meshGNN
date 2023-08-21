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

# Load Datasets

adjacencies, coords, areas, dataset_size = get_dataset('mesh800',True)
adjacencies_test, coords_test, areas_test, dataset_size_test = get_dataset('mesh200',True)


epochs=10
batch_sizes = [4,8,16,32,64]
epoch_loss_container = []
epoch_loss_val_container = []

# Create training log file

for batch_size in batch_sizes:

  model = SageBase(64,32,3,2).to(device)

  optimizer = torch.optim.Adam(model.parameters(), lr=1e-4, weight_decay=1e-5)
  scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',factor=0.5, patience=20, threshold=0.0001, threshold_mode='abs')

  loss_array = []
  loss_val_array = []
  epoch_loss_array = []
  epoch_loss_val_array = []
  lr_array = []

  start = time.time()
  
  print(f'\nTraining started:  batch_size = {batch_size}\n')
  
  for epoch in range(1, epochs+1):

      adjacencies_s, coords_s, areas_s = shuffle(adjacencies, coords, areas)
      adjacencies_test_s, coords_test_s, areas_test_s = shuffle(adjacencies_test, coords_test, areas_test)
    
      model.train()
      loss=torch.tensor([0.]).to(device)

      # Training
      for i in range(dataset_size):
        g = (get_sample(adjacencies_s, coords_s, areas_s,i,
                        randomRotate=bool(random.getrandbits(1)),selfloop=True)).to(device)
        out = model(g.x,g.edge_index) 
        loss += loss_normalized_cut(out, g)
        
        if i%batch_size==0 or i==dataset_size-1:
          optimizer.zero_grad()  
          loss_array.append(loss.item())
          loss.backward()  
          optimizer.step()  
          loss=torch.tensor([0.]).to(device)
      
      epoch_loss_array.append(sum(loss_array[-batch_size:])/batch_size)
      scheduler.step(loss_array[-1])
      lr_array.append(optimizer.param_groups[0]['lr'])

      # Validation  
      model.eval()
      loss_val = torch.tensor([0.]).to(device)

      for i in range(dataset_size_test):
        g_val = (get_sample(adjacencies_test, coords_test, areas_test,i,
                            randomRotate=bool(random.getrandbits(1)),selfloop=True)).to(device)
        out_val = model(g_val.x,g_val.edge_index) 
        loss_val += loss_normalized_cut(out_val, g_val)
        if i%batch_size==0 or i==dataset_size_test-1:
          loss_val_array.append(loss_val.item())
          loss_val = torch.tensor([0.]).to(device)
      epoch_loss_val_array.append(sum(loss_val_array[-batch_size:])/batch_size)

      print(f'epoch: {epoch} \t\tloss: {epoch_loss_array[-1]:.5f}  \t\tval loss: {epoch_loss_val_array[-1]:.5f} \t\tlr: {lr_array[-1]}')
  
  epoch_loss_container.append(epoch_loss_array)
  epoch_loss_val_container.append(epoch_loss_val_array)

  training_time = (time.time()-start)/60

  print(f'\nTraining ended - runtime: {training_time:.2f} min\n')

np.savetxt("data/loss_batchsize.csv", np.asarray(epoch_loss_container), delimiter=",")
np.savetxt("data/loss_val_batchsize.csv", np.asarray(epoch_loss_val_container), delimiter=",")

# Plots (saved in outputs/plots)
# plot loss
plt.figure(figsize=(10,5))
for i in range(len(batch_sizes)):
  plt.plot(np.array(epoch_loss_container[i])/epoch_loss_container[i][0])  

plt.title('Hyperparameters search - batch size')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.legend(batch_sizes) 
plt.savefig("plots/loss_batchsize.png",dpi=300)

plt.figure(figsize=(10,5))
for i in range(len(batch_sizes)):
  plt.plot(np.array(epoch_loss_val_container[i])/epoch_loss_val_container[i][0])  

plt.title('Hyperparams search - Batch size')
plt.xlabel('Epochs')
plt.ylabel('Validation Loss')
plt.legend(batch_sizes) 
plt.savefig("plots/loss_val_batchsize.png",dpi=300)