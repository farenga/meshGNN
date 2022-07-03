%--------------------------------------------------------------------
% PURPOSE:
%
% This routine assembles the SIPG matrix.
%
% Author:
% Paola Antonietti & Marco Sarti
%--------------------------------------------------------------------

function [Matrices] = matrix2D(femregion,neighbour,Dati)

%GP: nodes and weights on the reference triangular element ((0,0) (1,0) (0,1))
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Dati.nqn);
nqn_1D = length(w_1D);

if femregion.fem == 0 % scaling the  penalty parameter as a function of the degree of the shape functions
    penalty_coeff = Dati.penalty_coeff;
else
    penalty_coeff = Dati.penalty_coeff.*(femregion.fem.^2); %GP: that is 10 * p^2;
end


V=sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (grad(u) grad(v) dx}
M=sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (u . v ) dx}
I=sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} {grad v} . [u]ds
S=sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} penalty  h_e^(-1) [v].[u] ds
f=sparse(femregion.ndof,1);               % \int_{\Omega} f . v dx + boundary conditions
index_shift=0;

for ie = 1:femregion.ne % loop over elements

    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift = index_element(end);

    neigh_ie = neighbour.neigh{ie};      %GP: list of the neigh element of ie, neigh_ie(i)=idx means idx is the neigh elem that share the edge i of ie
    neighedges_ie = neighbour.neighedges{ie}; %GP: neighedges_ie(i)=j means that j is the edge of idx shared with ie, whose edge is i
    coords_elem = femregion.coords_element{ie}; %GP: matrix 2xNvertex of coordinate verteces of ie
    
    [normals,meshsize] = get_normals_meshsize_faces(coords_elem);
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    for iTria = 1:size(Tria,1)

        v1 = coords_elem(Tria(iTria,1),:);
        v2 = coords_elem(Tria(iTria,2),:);
        v3 = coords_elem(Tria(iTria,3),:);
        
        %GP: pphys_2D will be the nodes on the triangular (v1 v2 v3), the
        %map BJ is from ((0,0) (1,0) (0,1)) to ((x0,y0) (x1,y1) (x2,y2))
        [BJ, BJinv, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        Jdet=det(BJ);                       % determinant

        [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);

        for k = 1:length(w_2D) % loop over 2D quadrature nodes

            dx = w_2D(k)*Jdet;

            x = pphys_2D(k,1);
            y = pphys_2D(k,2);
            F = eval(Dati.source);
            for i = 1:femregion.nln % loop over scalar shape functions
                f(index(i)) = f(index(i))+F*dphiq(k,i).*dx;
                for j = i:femregion.nln % loop over scalar shape functions
                    V(index(i),index(j)) = V(index(i),index(j)) + Grad(k,:,i) * Grad(k,:,j)' .*dx ;
                    M(index(i),index(j)) = M(index(i),index(j)) + dphiq(k,i) * dphiq(k,j) .*dx ;
%                     if i~=j
%                         M(index(j),index(i)) = M(index(j),index(i)) + dphiq(k,i) * dphiq(k,j) .*dx ;
%                         V(index(j),index(i)) = V(index(j),index(i)) + Grad(k,:,i) * Grad(k,:,j)' .*dx ;
%                     end
                      if k==length(w_2D)
                        M(index(j),index(i)) = M(index(i),index(j));
                        V(index(j),index(i)) = V(index(i),index(j));
                    end
                end
            end
        end
    end
    
    IN = zeros(femregion.nln,femregion.nln,neighbour.nedges(ie));
    SN = zeros(femregion.nln,femregion.nln,neighbour.nedges(ie));
    for iedg = 1:neighbour.nedges(ie) % loop over faces, neighbour.nedges(ie) is the number of edges of the element ie

        neigedge = neighedges_ie(iedg);    % index of neighbour edge
 
        %%%%%%%%%%% Scaling: paper Cangiani, Georgoulis, Houston
        Cinv = femregion.area(ie)./femregion.max_kb{ie};
        
        if neigh_ie(iedg) == -1 %GP: means that iedg is a boundary edge of ie
            penalty_scaled = penalty_coeff*Cinv(iedg)*meshsize(iedg)/femregion.area(ie);
        else
            Cinv_ext = femregion.area(neigh_ie(iedg))./femregion.max_kb{neigh_ie(iedg)}(neigedge);
            s1 = penalty_coeff*Cinv(iedg)*meshsize(iedg)/femregion.area(ie);
            s2 = penalty_coeff*Cinv_ext*meshsize(iedg)/femregion.area(neigh_ie(iedg));
            penalty_scaled = max([s1 s2]);
        end
        %%%%%%%%%%%%%%%%%%%%

        if iedg < neighbour.nedges(ie)
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(iedg+1,:)';
            mean = 0.5*(p1+p2);
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        else
            p1 = coords_elem(iedg,:)'; p2 = coords_elem(1,:)';
            mean = 0.5*(p1+p2);
            vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy]; %GP: v = [0,0] ??
            v_hat = [-vy;vx];
            p3 = mean+v_hat;
            v = [p1';p2';p3'];
        end            
    
        [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);
        
        [B_edge,G_edge] = evalshape2D(femregion, ie ,pphys_1D);
        
        if neigh_ie(iedg) ~= -1
              [B_edge_neigh] = evalshape2D(femregion, neigh_ie(iedg),pphys_1D);
        end
        
        for k = 1:nqn_1D   % loop over 1D quadrature nodes

            ds = meshsize(iedg)*w_1D(k);
            
            for i = 1:femregion.nln % loop over scalar shape functions

                for j = 1:femregion.nln % loop over scalar shape functions

                    S(index(i),index(j)) = S(index(i),index(j)) + penalty_scaled .* B_edge(k,i) .* B_edge(k,j) .* ds;
    
                    if neigh_ie(iedg) ~= -1 % internal faces

                        I(index(i),index(j)) = I(index(i),index(j))  +  0.5 .* ((G_edge(k,:,i))*normals(:,iedg)) .* B_edge(k,j) .* ds;
                        IN(i,j,iedg)		 = IN(i,j,iedg)          -  0.5 .* ((G_edge(k,:,i))*normals(:,iedg)) .* B_edge_neigh(k,j) .* ds;
                        SN(i,j,iedg)         = SN(i,j,iedg)          - penalty_scaled .* B_edge(k,i) .* B_edge_neigh(k,j) .* ds;

                    elseif neigh_ie(iedg) == -1 % boundary faces
                        I(index(i),index(j))  = I(index(i),index(j)) +  ((G_edge(k,:,i))*normals(:,iedg)) .* B_edge(k,j) .* ds;
                    end
                end

                if  neigh_ie(iedg) == -1 % boundary conditions
                    x = pphys_1D(k,1);
                    y = pphys_1D(k,2);
                    gd = eval(Dati.exact_sol);
                    f(index(i)) = f(index(i))  + penalty_scaled .* B_edge(k,i) .* gd .* ds ;
                    f(index(i)) = f(index(i))  - ((G_edge(k,:,i))*normals(:,iedg)) .* gd  .* ds;
                end
            end
        end
    end
    [I] = assemble_neigh(I,index,neigh_ie,IN,femregion.nln,neighbour.nedges(ie)) ; % assemble the neighbours local matrices
    [S] = assemble_neigh(S,index,neigh_ie,SN,femregion.nln,neighbour.nedges(ie)); % assemble the neighbours local matrices
end

Matrices = struct('A',V -transpose(I) -I +S,...
    'V',V,...
    'M',M,...
    'I',I,...
    'S',S,...
    'f',f);



