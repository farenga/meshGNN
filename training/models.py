import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

class SageBase(nn.Module):
    def __init__(self, hidden_units, lin_hidden_units, num_features, out_classes):
        super().__init__()
        
        self.conv1 = SAGEConv(num_features, hidden_units,aggr='mean')  
        self.conv2 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        self.conv3 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        self.conv4 = SAGEConv(hidden_units, hidden_units,aggr='mean')  
        
        self.lin1 = nn.Linear(hidden_units, lin_hidden_units)
        self.lin2 = nn.Linear(lin_hidden_units, lin_hidden_units)
        self.lin_last = nn.Linear(lin_hidden_units, out_classes)
        
        self.act = torch.tanh

    def normalize(self,x):
        coords_sample = x[:,:2]
        areas_sample = x[:,-1].unsqueeze(-1)

        max_coords = (torch.max(coords_sample,0)).values
        min_coords = (torch.min(coords_sample,0)).values
        
        if((max_coords[1]-min_coords[1])>(max_coords[0]-min_coords[0])):
            theta = torch.tensor(torch.pi/2).to(device)
            R = torch.tensor([[torch.cos(theta),-torch.sin(theta)],[torch.sin(theta),torch.cos(theta)]]).to(device)
            coords_sample = torch.matmul(R,coords_sample.t()).t()
        
        coords_sample = (coords_sample-torch.mean(coords_sample,0))/(torch.max(coords_sample,0)).values
        areas_sample = areas_sample/torch.max(areas_sample,0).values

        return torch.cat([coords_sample,areas_sample],-1)
      
    def forward(self, x, edge_index):
        
        x = self.normalize(x)

        x = self.act(self.conv1(x, edge_index))
        x = self.act(self.conv2(x, edge_index))
        x = self.act(self.conv3(x, edge_index))
        x = self.act(self.conv4(x, edge_index))

        x = self.act(self.lin1(x))
        x = self.act(self.lin2(x))
        x = self.lin_last(x)

        x = F.softmax(x,dim=1)

        return x


class SageRes(nn.Module):
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
            theta = torch.tensor(torch.pi/2).to(device)
            R = torch.tensor([[torch.cos(theta),-torch.sin(theta)],[torch.sin(theta),torch.cos(theta)]]).to(device)
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