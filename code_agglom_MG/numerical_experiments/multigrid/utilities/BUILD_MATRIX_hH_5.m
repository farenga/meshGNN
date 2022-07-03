function M_hH = BUILD_MATRIX_hH_5(femregion_H,femregion_h,nqn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% BUILD MATRIX M_hH %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_hH = sparse(femregion_h.ndof,femregion_H.ndof);

%Quadrature nodes on the referent domain
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(nqn);
nqn_1D = length(w_1D);

%[connect_edges, total_edges] = CREATE_CONNECTIVITY_EDGES(femregion_h);
%[connect_nodes] = CREATE_CONNECTIVITY_NODES(femregion_h);

for ie = 1:femregion_H.ne
    
    %index in sparse matrix
    index_H = (ie-1)*femregion_H.nln*ones(femregion_H.nln,1) + [1:femregion_H.nln]';
    
    %S(1).P(1).x    = [femregion_H.coords_element{ie}(:,1);femregion_H.coords_element{ie}(1,1)];
    %S(1).P(1).y    = [femregion_H.coords_element{ie}(:,2);femregion_H.coords_element{ie}(1,2)];
    %S(1).P(1).hole = 0;
    
    P1.x    = [femregion_H.coords_element{ie}(:,1);femregion_H.coords_element{ie}(1,1)];
    P1.y    = [femregion_H.coords_element{ie}(:,2);femregion_H.coords_element{ie}(1,2)];
    P1.hole = 0;
    
    xH_min = femregion_H.BBox(ie,1); xH_max = femregion_H.BBox(ie,2);
    yH_min = femregion_H.BBox(ie,3); yH_max = femregion_H.BBox(ie,4);
    
    a = union( find( femregion_h.BBox(:,1) > xH_max ), find( femregion_h.BBox(:,2) < xH_min ) );
    b = union( find( femregion_h.BBox(:,3) > yH_max ), find( femregion_h.BBox(:,4) < yH_min ) );
    
    inter_element = setdiff( [1:1:femregion_h.ne], union(a,b) );
    clear a, clear b;
    
    %find intersections
    for il = inter_element
        
        index_h = (il-1)*femregion_h.nln*ones(femregion_h.nln,1) + [1:femregion_h.nln]';
        
        %S(2).P(1).x    = [femregion_h.coords_element{il}(:,1);femregion_h.coords_element{il}(1,1)];
        %S(2).P(1).y    = [femregion_h.coords_element{il}(:,2);femregion_h.coords_element{il}(1,2)];
        %S(2).P(1).hole = 0;
        
        P2.x    = [femregion_h.coords_element{il}(:,1);femregion_h.coords_element{il}(1,1)];
        P2.y    = [femregion_h.coords_element{il}(:,2);femregion_h.coords_element{il}(1,2)];
        P2.hole = 0;
        
        %tstart = tic;
        %P = PolygonClip(S(1).P,S(2).P,1);
        Polygon_intersection = PolygonClip(P1,P2,1);
        %telapse = toc(tstart);
        %X = ['Polyclip: ',num2str(telapse),'.'];
        %disp(X);
        
        if size(Polygon_intersection,2)~=0
            
            % When the intersection is composed by more than one domain
            for idx_intersec = 1:size(Polygon_intersection,2)
                P = Polygon_intersection(idx_intersec);
                nedges = length(P.x);
                
                
                %Build Delaunay triangulation of the Polygon P
                coords_elem = [P.x, P.y];
                edges = [[1:nedges]' [2:nedges 1]'];
                Tria_Del = DelaunayTri(P.x, P.y, edges);
                io = Tria_Del.inOutStatus();
                Tria = Tria_Del.Triangulation(io==1,:);
                
                
                %Make integrals
                for iTria = 1:size(Tria,1)
                    v1 = coords_elem(Tria(iTria,1),:);
                    v2 = coords_elem(Tria(iTria,2),:);
                    v3 = coords_elem(Tria(iTria,3),:);
                    
                    [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
                    Jdet=det(BJ);
                    
                    [dphiq_H,Grad_H] = evalshape2D(femregion_H, ie,pphys_2D);
                    [dphiq_h,Grad_h] = evalshape2D(femregion_h, il,pphys_2D);
                    
                    for k = 1:length(w_2D) % loop over 2D quadrature nodes
                        
                        dx = w_2D(k)*Jdet;
                        
                        for i = 1:femregion_h.nln % loop over scalar shape functions of V_h
                            for j = 1:femregion_H.nln % loop over scalar shape functions of V_H
                                M_hH(index_h(i),index_H(j)) = M_hH(index_h(i),index_H(j)) + dphiq_h(k,i) * dphiq_H(k,j) .*dx ;
                            end
                        end
                    end
                    
                end
            end
            
        end
    end
    
    %Delete edge if it is fully inside the polygon
    %if length(index_to_delete) ~= 0
    %    connect_edges(:,index_to_delete) = [];
    %end
end



end %end function