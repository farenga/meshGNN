function [M] = QuadFreeOnlyMass2D(femregion,neighbour,P,dP,Dati)

p = Dati.fem;
toll = 0.1; % this is the tololerance used by the IntegrationOverPolygon function, it needs to be choosen accurately!

% \Omega is the domain, E_h is the set containing the faces of the mesh
%V = sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (grad(u) grad(v) dx}
M = sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (u . v ) dx}
%I = sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} {grad v} . [u]ds
%S = sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} penalty  h_e^(-1) [v].[u] ds
%f = sparse(femregion.ndof,1);               % \int_{\Omega} f . v dx + boundary conditions



% OFF LINE COMPUTATIONS:

% 1) local enumeration basis functions:
%    - Phi_1 = Phi_(0,0)
%    - Phi_2 = Phi_(0,1)
%    - Phi_3 = Phi_(0,2)
%    ...
%    - Phi_{p+1} = Phi_(0,p)
%    - Phi_{p+2} = Phi_(1,0)
%    ...
%    - Phi_{2*p+1} = Phi_(1,p-1)
%    ...
%    - Phi_{femregion.nln} = Phi_(p,0)
local_enumeration_basis_functions = zeros(femregion.nln,2);
idx = 1;
for ix = 0:p
    for iy = 0:p-ix
        local_enumeration_basis_functions(idx,1) = ix;
        local_enumeration_basis_functions(idx,2) = iy;
        idx = idx + 1;
    end
end


%index_shift = 0;



% 3) Partial sum of coefficient in reference element:
%     this routine receive P and dP, coefficient matrices such that:
%     phi^_{ix}(x)  =  \sum_{i=0,...,ix} P(ix,i)  * x^i
%     dphi^_{ix}(x) =  \sum_{i=0,...,ix} dP(ix,i) * x^i
%     that is the homogeneous expansion of the 1D reference shape function
%     and its derivative
%
%    Partial sum of shape functions' coefficients in the reference element
%     - C_x1_x2(ix,jx,k)   =  \sum_{i+j=k} P(ix,i) * P(jx,j)
%     - C_dx1_x2(ix,jx,k)  =  \sum_{i+j=k} dP(ix,i) * P(jx,j)
%     - C_x1_dx2(ix,jx,k)  =  \sum_{i+j=k} P(ix,i) * dP(jx,j)
%     - C_dx1_dx2(ix,jx,k) =  \sum_{i+j=k} dP(ix,i) * dP(jx,j)
C_x1_x2 = zeros(p+1,p+1,2*p+1);  C_dx1_dx2 = zeros(p+1,p+1,2*p+1);
C_x1_dx2 = zeros(p+1,p+1,2*p+1);  C_dx1_x2 = zeros(p+1,p+1,2*p+1);
for idx = 0:p
    for idy = idx:p
        
        for k = 0:idx+idy
            Cxy = 0; Cdxdy = 0; Cdxy = 0; Cxdy = 0;
            for i=0:min(k,idx)
                if 0 <= k-i && k-i <= idy
                    Cxy = Cxy + P(idx+1,i+1)*P(idy+1,k-i+1);
                    Cdxdy = Cdxdy + dP(idx+1,i+1)*dP(idy+1,k-i+1);
                    Cdxy = Cdxy + dP(idx+1,i+1)*P(idy+1,k-i+1);
                    Cxdy = Cxdy + P(idx+1,i+1)*dP(idy+1,k-i+1);
                end
            end
            C_x1_x2(idx+1,idy+1,k+1) = Cxy;
            C_dx1_dx2(idx+1,idy+1,k+1) = Cdxdy;
            C_dx1_x2(idx+1,idy+1,k+1) = Cdxy;
            C_x1_dx2(idx+1,idy+1,k+1) = Cxdy;
            
            % thanks to some symmetric properties
            if idy ~= idx
                C_x1_x2(idy+1,idx+1,k+1) = Cxy;
                C_dx1_dx2(idy+1,idx+1,k+1) = Cdxdy;
                C_dx1_x2(idy+1,idx+1,k+1) = Cxdy;
                C_x1_dx2(idy+1,idx+1,k+1) = Cdxy;
            end
        end
        
    end
end


for ie = 1:femregion.ne % loop over elements
    
    % with the index ie, we are selecting the general element \kappa of the
    % mesh
    
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    %neigh_ie = neighbour.neigh{ie};      % list of the neigh element of ie, neigh_ie(i)=idx means idx is the neigh elem that share the edge i of ie
    %neighedges_ie = neighbour.neighedges{ie}; % neighedges_ie(i)=j means that j is the edge of idx shared with ie, whose edge is i
    v = femregion.coords_element{ie}';
    
    %[normals,meshsize] = get_normals_meshsize_faces(femregion.coords_element{ie});
    
    BBox = femregion.BBox(ie,:);
    x1B = BBox(1); x2B = BBox(2);
    y1B = BBox(3); y2B = BBox(4);
    
    x0 = x1B;   % x-coordinates of BBox vertices
    x1 = x2B;
    x2 = x2B;
    x3 = x1B;
    
    y0 = y1B;   % y-coordinates of BBox vertices
    y1 = y1B;
    y2 = y2B;
    y3 = y2B;
    
    % define Map F(x^,y^) = BJ*[x^,y^]' + trans
    trans = (0.25) .* [ x0 + x1 + x2 + x3 ;  y0 + y1 + y2 + y3];
    BJ = .25 .* [-x0 + x1 + x2 - x3 ,  -x0 - x1 + x2 + x3 ; -y0 + y1 + y2 - y3 , -y0 - y1 + y2 + y3];
    Jdet = det(BJ);
    
    % inverse of F
    D = BJ(2,2)*BJ(1,1)-BJ(1,2)*BJ(2,1);
    BJ_inv = [BJ(2,2) -BJ(1,2); -BJ(2,1) BJ(1,1)]/D;
    trans_inv = [-BJ(2,2)*trans(1)+BJ(1,2)*trans(2); BJ(2,1)*trans(1)-BJ(1,1)*trans(2)]/D;
    
    %coeff_a = BJ_inv(1,1)*BJ_inv(1,1);
    %coeff_b = BJ_inv(2,2)*BJ_inv(2,2);
    
    % Define \kappa^ = Finv(\kappa)
    vhat = [];
    for i = 1:size(v,2)
        vhat = [vhat, BJ_inv*v(:,i) + trans_inv];
    end
    
    % Homoenenous mass matrix on \kappa^ = F_inv(\kappa):
    % HM(k+1,l+1) = \int_{\kappa^} x^k * y^l dxdy;
    % FHM(k+1,l+1, nface) = \int_{E_nface^} x^k * y^l; E_nface^ is the
    % nface-th of \kappa^
    HM = zeros(p+1,p+1); FHM = zeros(p+1,p+1,neighbour.nedges(ie));
    for k = 0:2*p
        for l = 0:2*p-k
            [HM(k+1,l+1), FHM(k+1,l+1,1:end)] = IntegralOverPolygonAndFaces(1, k, l, vhat, toll);
        end
    end
    
    % loop over shape functions
    for idI = 1:femregion.nln
        
        % Phi_I(x,y) = phi_{ix}(x) * phi_{iy}(y)
        ix = local_enumeration_basis_functions(idI,1);
        iy = local_enumeration_basis_functions(idI,2);
        
        for idJ = idI:femregion.nln % and then consider the symmetry of M and V
            
            jx = local_enumeration_basis_functions(idJ,1);
            jy = local_enumeration_basis_functions(idJ,2);
            
            I_mass = 0;
            
            %loop on homogeneous polynomial functions
            for k = 0:ix+jx
                for l = 0:iy+jy
                    
                    I_mass = I_mass + ...
                        C_x1_x2(ix+1,jx+1,k+1) * C_x1_x2(iy+1,jy+1,l+1) * Jdet * HM(k+1,l+1);
                    
                end
            end
            
            M(index(idI),index(idJ)) = I_mass;
            %V(index(idI),index(idJ)) = I_stif;
            
            % thanks to the symmetry of M and V
            M(index(idJ),index(idI)) = I_mass;
            
        end
    end
end


end % end function