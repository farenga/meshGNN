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


from models import SageBase, SageRes
from losses import loss_normalized_cut
from utils import get_dataset, get_sample

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Load Datasets

adjacencies, coords, areas, dataset_size = get_dataset('mesh800')
adjacencies_test, coords_test, areas_test, dataset_size_test = get_dataset('mesh200')

# Load model definition (from models.py)

#model = SageBase(64,32,3,2).to(device)
model = SageRes().to(device)

modelname = model.__class__.__name__
print('\n'+modelname)

print('# parameters: ',sum(p.numel() for p in model.parameters()))


# Train 

epochs=300
batch_size=4
optimizer = torch.optim.Adam(model.parameters(), lr=1e-5, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',factor=0.5, patience=20, threshold=0.0001, threshold_mode='abs')


loss_array = []
loss_val_array = []
epoch_loss_array = []
epoch_loss_val_array = []
lr_array = []


# Create training log file

logfile = Path('outputs/train_info/'+modelname+'_logfile.txt')
logfile.touch(exist_ok=True)
f = open(logfile,'w')

print(modelname+'\n',file=f)
start = time.time()

print('\nTraining started\n')
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
    print(f'epoch: {epoch} \t\tloss: {epoch_loss_array[-1]:.5f}  \t\tval loss: {epoch_loss_val_array[-1]:.5f} \t\tlr: {lr_array[-1]}',file=f)

training_time = (time.time()-start)/60

print(f'\nTraining ended - runtime: {training_time:.2f} min\n')
print(f'\nTraining runtime: {training_time:.2f} min\n',file=f)

# Save training infos (losses,lr schedule)
np.savetxt("outputs/train_info/loss_"+modelname+".csv", np.asarray(loss_array), delimiter=",")
np.savetxt("outputs/train_info/epoch_loss_"+modelname+".csv", np.asarray(epoch_loss_array), delimiter=",")
np.savetxt("outputs/train_info/loss_val_"+modelname+".csv", np.asarray(loss_val_array), delimiter=",")
np.savetxt("outputs/train_info/epoch_loss_val_"+modelname+".csv", np.asarray(epoch_loss_val_array), delimiter=",")
np.savetxt("outputs/train_info/lr_"+modelname+".csv", np.asarray(lr_array), delimiter=",")

# Save trained model
torch.save(model.state_dict(),"outputs/models/"+modelname+".pt")


# Plots (saved in outputs/plots)
# plot loss
plt.figure(figsize=(10,5))   
plt.plot(epoch_loss_array)  
plt.plot(epoch_loss_val_array) 
plt.xlabel('Epochs')
plt.legend(['Training loss', 'Validation loss']) 
plt.savefig("outputs/plots/loss"+modelname+".png",dpi=300)

# plot lr schedule
plt.figure(figsize=(10,5))   
plt.plot(lr_array)  
plt.xlabel('Epochs')
plt.ylabel('lr')
plt.savefig("outputs/plots/lr"+modelname+".png",dpi=300)
