function [Matrices] = QuadFreeMatrix2D_2(femregion,neighbour,P,dP,Dati)

p = Dati.fem;
toll = 0.1; % this is the tololerance used by the IntegrationOverPolygon function, it needs to be choosen accurately!

% \Omega is the domain, E_h is the set containing the faces of the mesh
V = sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (grad(u) grad(v) dx}
M = sparse(femregion.ndof,femregion.ndof);  % \int_{\Omega} (u . v ) dx}
I = sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} {grad v} . [u]ds
S = sparse(femregion.ndof,femregion.ndof);  % \int_{E_h} penalty  h_e^(-1) [v].[u] ds
f = sparse(femregion.ndof,1);               % \int_{\Omega} f . v dx + boundary conditions

% N.B.: this function admints only constant forcing term!

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

% 2) compute the penalty of the SIPG method
if femregion.fem == 0 % scaling the  penalty parameter as a function of the degree of the shape functions
    penalty_coeff = Dati.penalty_coeff;
else
    penalty_coeff = Dati.penalty_coeff.*(femregion.fem.^2); %GP: that is 10 * p^2;
end

index_shift = 0;



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


% 4) matrix of Newton's binomials
for i = 0:p
    for j = 0:i
        bin_coeff(i+1,j+1) = nchoosek(i,j);
    end
end


HM = zeros(2*p+1,2*p+1); FHM = zeros(2*p+1,2*p+1);
for ie = 1:femregion.ne % loop over elements
    
    % with the index ie, we are selecting the general element \kappa of the
    % mesh
    
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    neigh_ie = neighbour.neigh{ie};      % list of the neigh element of ie, neigh_ie(i)=idx means idx is the neigh elem that share the edge i of ie
    neighedges_ie = neighbour.neighedges{ie}; % neighedges_ie(i)=j means that j is the edge of idx shared with ie, whose edge is i
    v = femregion.coords_element{ie}';
    
    [normals,meshsize] = get_normals_meshsize_faces(femregion.coords_element{ie});
    
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
    
    coeff_a = BJ_inv(1,1)*BJ_inv(1,1);
    coeff_b = BJ_inv(2,2)*BJ_inv(2,2);
    
    % Define \kappa^ = Finv(\kappa)
    vhat = [];
    for i = 1:size(v,2)
        vhat = [vhat, BJ_inv*v(:,i) + trans_inv];
    end
    
    % Integral over interfaces
    IN = zeros(femregion.nln,femregion.nln,neighbour.nedges(ie));
    SN = zeros(femregion.nln,femregion.nln,neighbour.nedges(ie));
    
    % Homoenenous mass (HM) integration on \kappa^ = F_inv(\kappa): and on faces
    % HM(k+1,l+1) = \int_{\kappa^} x^k * y^l dxdy;
    % FHM(k+1,l+1, nface) = \int_{E_nface^} x^k * y^l; E_nface^ is the
    % nface-th of \kappa^
    HM = 0*HM;
    d = size(vhat,1);
    m = size(vhat,2); % this is also equal to neighbour.nedges(ie)
    for iedg = 1:neighbour.nedges(ie) % loop over faces, neighbour.nedges(ie) is the number of edges of the element ie
        
        % Here face integration, FHM is for Face Homoneneous functions' Mass matrix
        % FHM(k+1,l+1) = \int_{E_{iedg}^} x^k * y^l; E_{iedg}^ is the iedg-th face of \kappa^
        FHM = 0*FHM;
        if iedg < m
            x1 = vhat(1,iedg); x2 = vhat(1,iedg+1);
            y1 = vhat(2,iedg); y2 = vhat(2,iedg+1);
        else
            x1 = vhat(1,iedg); x2 = vhat(1,1);
            y1 = vhat(2,iedg); y2 = vhat(2,1);
        end
        
        bi = ((y2-y1)*x1 + (x1-x2)*y1) / sqrt( (y2-y1)^2 + (x1-x2)^2 );
        dij = sqrt( (y2-y1)^2 + (x2-x1)^2 );
        
        FHM(1,1) = dij/(d-1);
        HM(1,1) = HM(1,1) + bi*FHM(1,1)/d;
        
        for k=1:2*p
            FHM(k+1,1) = ( dij*x2^k + x1*k*FHM(k,1) )/(d+k-1);
            FHM(1,k+1) = ( dij*y2^k + y1*k*FHM(1,k) )/(d+k-1);
            
            HM(k+1,1) = HM(k+1,1) + bi*FHM(k+1,1)/(d+k);
            HM(1,k+1) = HM(1,k+1) + bi*FHM(1,k+1)/(d+k);
        end
        
        for k=1:2*p
            for l=1:2*p-k
                FHM(k+1,l+1) = ( dij*(x2^k)*(y2^l) + x1*k*FHM(k,l+1) + y1*l*FHM(k+1,l) )/(d+k+l-1);
                HM(k+1,l+1) = HM(k+1,l+1) + bi*FHM(k+1,l+1)/(d+k+l);
            end
        end
        % end of face integration on E_{iedg}^
        
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
        
        
        % Scaling for integration on the reference element
        if iedg < neighbour.nedges(ie)
            x1 = vhat(1,iedg); x2 = vhat(1,iedg+1);
            y1 = vhat(2,iedg); y2 = vhat(2,iedg+1);
        else
            x1 = vhat(1,iedg); x2 = vhat(1,1);
            y1 = vhat(2,iedg); y2 = vhat(2,1);
        end
        
        if x1 ~= x2
            mhat = (y2 - y1) / (x2 - x1);
            JdetFace = sqrt( (BJ(1,1)^2 + (mhat*BJ(2,2))^2) / (1+mhat^2) );
        else
            JdetFace = BJ(2,2);
        end
        
        
        % from here ap to the line 425, there are the instruction to evaluate
        % the integral 
        %    \int_{E_h} penalty  h_e^(-1) [Phi_I].[Phi_J] ds
        % when supp(Phi_I) = \kappa, and supp(Phi_J)= \kappa' \ne \kappa
        Btilde11 = 0; Btilde22 = 0; ttilde1 = 0; ttilde2 = 0;
        Ctilde1 = zeros(p+1,p+1); Ctilde2 = zeros(p+1,p+1);
        Ctilde1_x1_x2 = zeros(p+1,p+1,2*p+1); Ctilde2_x1_x2 = zeros(p+1,p+1,2*p+1);
        Ctilde1_dx1_x2 = zeros(p+1,p+1,2*p+1); Ctilde2_dx1_x2 = zeros(p+1,p+1,2*p+1);
        
        if neigh_ie(iedg) ~= -1
            
            % Map of the neighnour
            BBox_ne = femregion.BBox(neigh_ie(iedg),:);
            x1B = BBox_ne(1); x2B = BBox_ne(2);
            y1B = BBox_ne(3); y2B = BBox_ne(4);
            
            x0 = x1B;   % x-coordinates of BBox vertices
            x1 = x2B;
            x2 = x2B;
            x3 = x1B;
            
            y0 = y1B;   % y-coordinates of BBox vertices
            y1 = y1B;
            y2 = y2B;
            y3 = y2B;
            
            % define Map F(x^,y^) = BJ*[x^,y^]' + trans
            trans_ne = (0.25) .* [ x0 + x1 + x2 + x3 ;  y0 + y1 + y2 + y3];
            BJ_ne = .25 .* [-x0 + x1 + x2 - x3 ,  -x0 - x1 + x2 + x3 ; -y0 + y1 + y2 - y3 , -y0 - y1 + y2 + y3];
            Jdet_ne = det(BJ_ne);
            
            % inverse of F
            %D = BJ(2,2)*BJ(1,1)-BJ(1,2)*BJ(2,1);
            BJ_inv_ne = [BJ_ne(2,2) -BJ_ne(1,2); -BJ_ne(2,1) BJ_ne(1,1)]/Jdet_ne;
            trans_inv_ne = [-BJ_ne(2,2)*trans_ne(1) + BJ_ne(1,2)*trans_ne(2);...
                BJ_ne(2,1)*trans_ne(1) - BJ_ne(1,1)*trans_ne(2)]/Jdet_ne;
            
            Btilde11 = BJ_inv_ne(1,1)*BJ(1,1); Btilde22 = BJ_inv_ne(2,2)*BJ(2,2);
            ttilde1 = BJ_inv_ne(1,1)*trans(1) + trans_inv_ne(1);
            ttilde2 = BJ_inv_ne(2,2)*trans(2) + trans_inv_ne(2);
            
            
            % Local coeff for neighbour shape functions
            for k = 0:p
                for j = k:p
                    
                    for l = k:j
                        Ctilde1(j+1,k+1) = Ctilde1(j+1,k+1) + ...
                            P(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde11^k * ttilde1^(l-k);
                        
                        Ctilde2(j+1,k+1) = Ctilde2(j+1,k+1) + ...
                            P(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde22^k * ttilde2^(l-k);
                    end
                    
                end
            end
            
            
            for idx = 0:p
                for idy = 0:p
                    
                    for k = 0:idx+idy
                        Cxy1 = 0; Cxy2 = 0; Cdxy1 = 0; Cdxy2 = 0;
                        for i = max(0,k-idy):min(k,idx)
                            Cxy1 = Cxy1 + P(idx+1,i+1)*Ctilde1(idy+1,k-i+1);
                            Cxy2 = Cxy2 + P(idx+1,i+1)*Ctilde2(idy+1,k-i+1);
                            Cdxy1 = Cdxy1 + dP(idx+1,i+1)*Ctilde1(idy+1,k-i+1);
                            Cdxy2 = Cdxy2 + dP(idx+1,i+1)*Ctilde2(idy+1,k-i+1);
                        end
                        Ctilde1_x1_x2(idx+1,idy+1,k+1) = Cxy1;
                        Ctilde2_x1_x2(idx+1,idy+1,k+1) = Cxy2;
                        Ctilde1_dx1_x2(idx+1,idy+1,k+1) = Cdxy1;
                        Ctilde2_dx1_x2(idx+1,idy+1,k+1) = Cdxy2;
                    end
                    
                end
            end
            
        end
        
        % loop over shape functions
        for idI = 1:femregion.nln
            
            ix = local_enumeration_basis_functions(idI,1);
            iy = local_enumeration_basis_functions(idI,2);
            
            for idJ = 1:femregion.nln
                
                jx = local_enumeration_basis_functions(idJ,1);
                jy = local_enumeration_basis_functions(idJ,2);
                
                Int_S = 0; Int_I = 0;
                Int_SN = 0; Int_IN = 0;
                
                %loop on homogeneous polynomial functions
                for k = 0:ix+jx
                    for l = 0:iy+jy
                        
                        Int_S = Int_S + C_x1_x2(ix+1,jx+1,k+1) * C_x1_x2(iy+1,jy+1,l+1) ...
                            * JdetFace * FHM(k+1,l+1);
                        
                        Ca = C_dx1_x2(ix+1,jx+1,k+1)*C_x1_x2(iy+1,jy+1,l+1)*normals(1,iedg)*BJ_inv(1,1);
                        Cb = C_x1_x2(ix+1,jx+1,k+1)*C_dx1_x2(iy+1,jy+1,l+1)*normals(2,iedg)*BJ_inv(2,2);
                        Int_I = Int_I + (Ca + Cb) * JdetFace * FHM(k+1,l+1);
                        
                    end
                end
                
                if neigh_ie(iedg) ~= -1
                    for k = 0:ix+jx
                        for l = 0:iy+jy
                            
                            Int_SN = Int_SN + Ctilde1_x1_x2(ix+1,jx+1,k+1) * Ctilde2_x1_x2(iy+1,jy+1,l+1) ...
                                * JdetFace * FHM(k+1,l+1);
                            
                            Ca = Ctilde1_dx1_x2(ix+1,jx+1,k+1)*Ctilde2_x1_x2(iy+1,jy+1,l+1)*normals(1,iedg)*BJ_inv(1,1);
                            Cb = Ctilde1_x1_x2(ix+1,jx+1,k+1)*Ctilde2_dx1_x2(iy+1,jy+1,l+1)*normals(2,iedg)*BJ_inv(2,2);
                            Int_IN = Int_IN + (Ca + Cb) * JdetFace * FHM(k+1,l+1);
                            
                        end
                    end
                end
                
                S(index(idI),index(idJ)) = S(index(idI),index(idJ)) + penalty_scaled * Int_S;
                
                if neigh_ie(iedg) ~= -1 % internal faces
                    I(index(idI),index(idJ)) = I(index(idI),index(idJ)) +  0.5 * Int_I;
                    IN(idI,idJ,iedg) = IN(idI,idJ,iedg) -  0.5 * Int_IN;
                    SN(idI,idJ,iedg) = SN(idI,idJ,iedg) - penalty_scaled * Int_SN;
                elseif neigh_ie(iedg) == -1 % boundary faces
                    I(index(idI),index(idJ)) = I(index(idI),index(idJ)) +  Int_I;
                end
                
            end
            
            if neigh_ie(iedg) == -1
                
                x = 0; y = 0; 
                g = eval(Dati.dirichlet);
                Int_S = 0;
                for k = 0:ix
                    for l = 0:iy
                        Int_S = Int_S + (penalty_scaled*P(ix+1,k+1)*P(iy+1,l+1) ...
                            - dP(ix+1,k+1)*P(iy+1,l+1)*normals(1,iedg) - P(ix+1,k+1)*dP(iy+1,l+1)*normals(2,iedg))*FHM(k+1,l+1);
                    end
                end
                f(index(idI)) = f(index(idI)) + g*JdetFace*Int_S;
            end
            
        end
        
    end
    I = assemble_neigh(I,index,neigh_ie,IN,femregion.nln,neighbour.nedges(ie)) ; % assemble the neighbours local matrices
    S = assemble_neigh(S,index,neigh_ie,SN,femregion.nln,neighbour.nedges(ie)); % assemble the neighbours local matrices
    
    % Again: this function only admints constant forcing term!
    x = 0; y = 0;
    ForcingHM = eval(Dati.source)*HM;
    
    % This is the computation of volume integrals
    % loop over shape functions
    for idI = 1:femregion.nln
        
        % Phi_I(x,y) = phi_{ix}(x) * phi_{iy}(y)
        ix = local_enumeration_basis_functions(idI,1);
        iy = local_enumeration_basis_functions(idI,2);
        
        for k = 0:ix
            for l = 0:iy
                f(index(idI)) = f(index(idI)) + P(ix+1,k+1)*P(iy+1,l+1)*ForcingHM(k+1,l+1)*Jdet;
            end
        end
        
        for idJ = idI:femregion.nln % and then consider the symmetry of M and V
            
            jx = local_enumeration_basis_functions(idJ,1);
            jy = local_enumeration_basis_functions(idJ,2);
            
            I_mass = 0; I_stif = 0;
            
            %loop on homogeneous polynomial functions
            for k = 0:ix+jx
                for l = 0:iy+jy
                    
                    I_mass = I_mass + ...
                        C_x1_x2(ix+1,jx+1,k+1) * C_x1_x2(iy+1,jy+1,l+1) * Jdet * HM(k+1,l+1);
                    
                    
                    Ca =  C_dx1_dx2(ix+1,jx+1,k+1) * C_x1_x2(iy+1,jy+1,l+1) * ...
                        coeff_a;
                    
                    Cb =  C_x1_x2(ix+1,jx+1,k+1) * C_dx1_dx2(iy+1,jy+1,l+1) * ...
                        coeff_b;
                    
                    I_stif = I_stif + (Ca + Cb) * Jdet * HM(k+1,l+1);
                    
                end
            end
            
            M(index(idI),index(idJ)) = I_mass;
            V(index(idI),index(idJ)) = I_stif;
            
            % thanks to the symmetry of M and V
            M(index(idJ),index(idI)) = I_mass;
            V(index(idJ),index(idI)) = I_stif;
            
        end
        
        
    end
    
end

Matrices = struct('A',V -transpose(I) -I +S,...
    'V',V,...
    'M',M,...
    'I',I,...
    'S',S,...
    'f',f);

end % end function