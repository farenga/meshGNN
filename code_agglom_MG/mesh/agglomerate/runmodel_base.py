import sys
import torch
import numpy as np
import torch_geometric
from torch_geometric.data import Data
import scipy.io
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
from torch_geometric.utils import softmax

class SageNet(torch.nn.Module):
    def __init__(self, hidden_units, lin_hidden_units, num_features, out_classes):
        super().__init__()
        
        self.conv1 = SAGEConv(num_features, hidden_units,aggr='mean')  
        self.conv2 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        self.conv3 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        self.conv4 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        
        self.lin1 = nn.Linear(hidden_units, lin_hidden_units)
        self.lin_last = nn.Linear(lin_hidden_units, out_classes)
        
        self.act = torch.tanh

    def normalize(self,x):
        coords_sample = x[:,:2]
        areas_sample = x[:,-1].unsqueeze(-1)

        max_coords = (torch.max(coords_sample,0)).values
        min_coords = (torch.min(coords_sample,0)).values
        
        if((max_coords[1]-min_coords[1])>(max_coords[0]-min_coords[0])):
            theta = torch.tensor(torch.pi/2)
            R = torch.tensor([[torch.cos(theta),-torch.sin(theta)],[torch.sin(theta),torch.cos(theta)]])
            coords_sample = torch.matmul(R,coords_sample.t()).t()
        
        coords_sample = (coords_sample-torch.mean(coords_sample,0))/(torch.max(coords_sample,0)).values
        areas_sample = areas_sample/torch.max(areas_sample,0).values

        return torch.cat([coords_sample,areas_sample],-1)
      
    def forward(self, x, edge_index):
        
        x = self.normalize(x)

        x = self.act(self.conv1(x, edge_index))
        x = self.act(self.conv2(x, edge_index))
        x = self.act(self.lin1(x))
        x = self.lin_last(x)
        x = F.softmax(x,dim=1)

        return x
        
        

def call_model(Adjacency, Coords, Areas):
    
    model = SageNet(64,32,3,2)
    model.load_state_dict(torch.load(sys.argv[1]))
    model.eval()
    
    A = torch.from_numpy(Adjacency)
    
    coords = torch.from_numpy(Coords)
    coords = coords.type(torch.float)
    areas = torch.from_numpy(Areas)
    areas = torch.unsqueeze(areas.type(torch.float),-1)
    x = torch.cat([coords,areas],-1)

    edge_index = (A > 0).nonzero().t()
    data = Data(x=x, edge_index=edge_index)
    
    out = model(data.x, data.edge_index)
    outclass = (out>=.5).int()
    
    return outclass.numpy()
    
z = call_model(Adjacency, Coords, Areas)