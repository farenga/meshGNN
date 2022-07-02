function IDX = aggl_metis_fun(mesh,weights,n_clust)

if nargin < 3
    n_clust = 2;
end

NE = mesh.elem_num(end); 
NN = mesh.elem_num(1);
eptr = ones(1,NE+1);
eind = [];
if nargin < 2
    weights = ones(1,NE);
end
for i = 1:NE
    elem = RView(mesh,i,mesh.dim);
    egs = edges(elem);
    vert_ids = elem.loc2glob{1}(unique(egs(:)'));
    eind = [eind,vert_ids];
    eptr(i+1) = eptr(i) + elem.elem_num(1);
    if nargin < 2
        polygon = edges2polygon(vertices(elem),edges(elem));
        P = polyshape(polygon);
        weights(i) = area(P);
    end
end


ncommon = 2; % 1 edge in common --> 2 nodes common in 2D
options.METIS_OPTION_DBGLVL = int64(2);  %
options.METIS_OPTION_UFACTOR = int64(100);  %
options = [];
% weights = [];
weights = int64(weights/min(weights)); % <--!!! weights must be interger (32 or 64 bits)



IDX = METIS_PartMeshDual(...
    int64(NE),int64(NN),int64(eptr),int64(full(eind)),...
    int64(weights),[],int64(ncommon),int64(n_clust),[],options);

if all(IDX==IDX(1))
    % try with no weights
    IDX = METIS_PartMeshDual(...
    int64(NE),int64(NN),int64(eptr),int64(full(eind)),...
    [],[],int64(ncommon),int64(n_clust),[],options);
end

if all(IDX==IDX(1))
    IDX = aggl_kmeans_fun(mesh,weights);
end

% if all(IDX==IDX(1))
%     error('failed partitioning with metis')
% end

end