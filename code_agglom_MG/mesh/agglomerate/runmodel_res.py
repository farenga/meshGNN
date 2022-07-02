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

class SageUNet(torch.nn.Module):
    def __init__(self):
        super().__init__()
        
        self.conv1 = SAGEConv(3, 8)  
        self.conv2 = SAGEConv(8, 16)
        self.conv3 = SAGEConv(16, 32)
        self.conv4 = SAGEConv(32, 64)
        
        self.conv1r = SAGEConv(64, 32)  
        self.conv2r = SAGEConv(32, 16)
        self.conv3r = SAGEConv(16, 8)
        
        self.lin1 = nn.Linear(8, 32)
        self.lin2 = nn.Linear(32, 16)
        self.lin_last = nn.Linear(16, 2)
        
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

        x1 = self.act(self.conv1(x, edge_index))
        x2 = self.act(self.conv2(x1, edge_index))
        x3 = self.act(self.conv3(x2, edge_index))
        x4 = self.act(self.conv4(x3, edge_index))

        x = self.act(self.conv1r(x4, edge_index))+x3
        x = self.act(self.conv2r(x, edge_index))+x2
        x = self.act(self.conv3r(x, edge_index))+x1

        x = self.act(self.lin1(x))
        x = self.act(self.lin2(x))
        x = self.lin_last(x)
    
        x = F.softmax(x,dim=1)

        return x
        

def call_model(Adjacency, Coords, Areas):
    
    model = SageUNet()
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