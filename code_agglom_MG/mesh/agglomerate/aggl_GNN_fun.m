function IDX = aggl_GNN_fun(mesh,weights,net)

% weights è un vettore che contiene l'area di ogni poligono
% IDX è un vettore che contine il cluster (1 o 2) di ogni poligono

% ATTENTION: Copy the directory to the model file
model_directory = 'C:\Users\gabri\Desktop\ProgettoNAPDE\meshGNN\code_agglom_MG\mesh\agglomerate\model_base.pt';


[graph,W] = connectivity(mesh);

Adjacency = W>0;
Areas = area_faces(mesh)';
Coords = cell2mat(graph.elem(1));

Adjacency = py.numpy.array(Adjacency);
Areas = py.numpy.array(Areas);
Coords = py.numpy.array(Coords);

y = pyrunfile(['runmodel_base.py ',model_directory], "z", Adjacency=Adjacency, Areas=Areas, Coords=Coords);
y = int64(y);
IDX = y(:,1) + 1;

% In case the model returns only one cluster instead of two
if all(IDX == IDX(1))
    IDX = aggl_kmeans_fun(mesh,weights);
end

end