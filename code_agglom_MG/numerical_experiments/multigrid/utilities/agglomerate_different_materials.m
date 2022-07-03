function [TauFinal] = agglomerate_different_materials(Tau_fine,Nmaterials)


for k = 1:Nmaterials
    
    % Get triangular local region related to the k-th material (green)
    I = find(Tau_fine.tag==k);
    Node = Tau_fine.coord;
    Element = {};
    for j=1:length(I)
        Element{j} = Tau_fine.connectivity{I(j)};
    end
    [Tau] = generate_local_region(Node,Element);
    %neigh = neighbours(Tau);
    %[error] = Plot_PolyMesh(Tau,1,'k',neigh,'-','b');
    a=1;
    
    
    neigh = neighbours(Tau);
    TauLoc{k} = agglomeration_metis(Tau,neigh,ceil(Tau.ne));
    
    %neigh = neighbours(Tau);
    %[error] = Plot_PolyMesh(Tau,1,'k',neigh,'-','b');
    
end

TauFinal = TauLoc{1}; % main mesh of material 1
for k = 2:Nmaterials
    
    starts = TauFinal.ne + 1;
    ends = TauFinal.ne + TauLoc{k}.ne;
    
    TauFinal.ne = TauFinal.ne + TauLoc{k}.ne;
    
    index = 1;
    for j=starts:ends
        TauFinal.connectivity{j} = TauLoc{k}.connectivity{index};
        TauFinal.coords_element{j} = TauLoc{k}.coords_element{index};
        TauFinal.edges_phys{j} = TauLoc{k}.edges_phys{index};
        TauFinal.edges{j} = TauLoc{k}.edges{index};
        Efinal(j) = k;
        index = index + 1;
    end
    
    TauFinal.nedges = [TauFinal.nedges; TauLoc{k}.nedges];
    TauFinal.centroids = [TauFinal.centroids; TauLoc{k}.centroids];
    TauFinal.BBox = [TauFinal.BBox; TauLoc{k}.BBox];
    TauFinal.edges_ID = [TauFinal.edges_ID; TauLoc{k}.edges_ID];
    a=1;
    
    
end

