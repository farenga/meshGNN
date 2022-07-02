import torch
import numpy as np
import scipy.io as sp
from torch_geometric.data import Data
from torch_geometric.utils import degree

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Data loading utilities

def get_dataset(dataset_name,relative=False):
  
  # load 
  if relative:
    adjacencies = np.squeeze(sp.loadmat('../datasets/'+dataset_name+'/AdjacencyMatrices.mat')['AdjacencyMatrices'],0)
    coords = np.squeeze(sp.loadmat('../datasets/'+dataset_name+'/CoordMatrices.mat')['CoordMatrices'],0)
    areas = np.squeeze(sp.loadmat('../datasets/'+dataset_name+'/AreaVectors.mat')['AreaVectors'],0)
  else:
    adjacencies = np.squeeze(sp.loadmat('datasets/'+dataset_name+'/AdjacencyMatrices.mat')['AdjacencyMatrices'],0)
    coords = np.squeeze(sp.loadmat('datasets/'+dataset_name+'/CoordMatrices.mat')['CoordMatrices'],0)
    areas = np.squeeze(sp.loadmat('datasets/'+dataset_name+'/AreaVectors.mat')['AreaVectors'],0)
    print('Dataset loaded')

  # clean
  rmv_idx = []
  for k in range(adjacencies.shape[-1]):
    if(np.sum(adjacencies[k])==0):
      rmv_idx.append(k)

  adjacencies = np.delete(adjacencies, rmv_idx)
  coords = np.delete(coords, rmv_idx)
  areas = np.delete(areas, rmv_idx)

  dataset_size = adjacencies.shape[-1]
  
  dataset_size = areas.shape[-1]
  if(rmv_idx):
      print('Dataset cleaned: removed', len(rmv_idx), 'null-arrays')

  print('Name:\t\t', dataset_name)
  print('Dimensions:\t',adjacencies.shape, areas.shape, coords.shape)


  return adjacencies, coords, areas, dataset_size

  
# returns a graph data structure sample from arrays of adj, coords, areas
def get_sample(adjacencies, coords, areas,i,randomRotate=False, selfloop=False, returnAll=False):

  coords_sample = torch.tensor(coords[i],dtype=torch.float)
  areas_sample = torch.tensor(areas[i],dtype=torch.float)

  if randomRotate:
    theta = torch.randn(1)*torch.pi
    #theta = torch.tensor(np.random.choice(np.array([torch.pi/2,0]),1),dtype=torch.float)
    R = torch.tensor([[torch.cos(theta),-torch.sin(theta)],[torch.sin(theta),torch.cos(theta)]])
    coords_rot = torch.matmul(R,coords_sample.t()).t()
    coords_sample = coords_rot

  x = torch.cat([coords_sample,areas_sample],-1)
  

  A = torch.tensor(adjacencies[i])+torch.eye(adjacencies[i].shape[0]) if selfloop else torch.tensor(adjacencies[i])
  edge_index = (A > 0).nonzero().t()

  data = Data(x=x, edge_index=edge_index)
  if returnAll:
    return data, A, x
  else:
    return data 


def rotate(coords):
    max_coords = (torch.max(coords,0)).values
    min_coords = (torch.min(coords,0)).values
    print(max_coords.shape)
    print(min_coords.shape)
    if((max_coords[1]-min_coords[1])>(max_coords[0]-min_coords[0])):
        theta = torch.tensor(torch.pi/2)
        R = torch.tensor([[torch.cos(theta),-torch.sin(theta)],[torch.sin(theta),torch.cos(theta)]])
        coords = torch.matmul(R,coords.t()).t()
    return coords
