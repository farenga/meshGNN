import torch
import numpy as np
import scipy.io as sp
from torch_geometric.data import Data
from torch_geometric.utils import degree

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# LOSS Definition
def loss_normalized_cut (y , graph):
  d = degree(graph.edge_index[0], num_nodes=y.size(0))
  gamma = torch.t(y) @ d
  c = torch.sum(y[graph.edge_index[0], 0]*y[graph.edge_index[1], 1])
  return torch.sum(torch.div(c,gamma)).to(device)