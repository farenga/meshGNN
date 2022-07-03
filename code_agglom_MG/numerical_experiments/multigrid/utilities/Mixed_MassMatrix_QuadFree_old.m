function M_hH = Mixed_MassMatrix_QuadFree(femregion_H,femregion_h,P_coeff,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% BUILD MATRIX M_hH %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%p = Dati.fem;
toll = 0.1;

M_hH = sparse(femregion_h.ndof,femregion_H.ndof);

local_enumeration_basis_functions = zeros(femregion_h.nln,2);
idx = 1;
for ix = 0:p
    for iy = 0:p-ix
        local_enumeration_basis_functions(idx,1) = ix;
        local_enumeration_basis_functions(idx,2) = iy;
        idx = idx + 1;
    end
end

for i = 0:p
    for j = 0:i
        bin_coeff(i+1,j+1) = nchoosek(i,j);
    end
end

for ie = 1:femregion_H.ne
    
    BBox = femregion_H.BBox(ie,:);
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
    
    %index in sparse matrix
    index_H = (ie-1)*femregion_H.nln*ones(femregion_H.nln,1) + [1:femregion_H.nln]';
    
    % polygon 1
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
        
        % polygon 2
        P2.x    = [femregion_h.coords_element{il}(:,1);femregion_h.coords_element{il}(1,1)];
        P2.y    = [femregion_h.coords_element{il}(:,2);femregion_h.coords_element{il}(1,2)];
        P2.hole = 0;
        
        %tstart = tic;
        Polygon_intersection = PolygonClip(P1,P2,1);
        
        Ctilde1 = zeros(p+1,p+1); Ctilde2 = zeros(p+1,p+1);
        Ctilde1_x1_x2 = zeros(p+1,p+1,2*p+1); Ctilde2_x1_x2 = zeros(p+1,p+1,2*p+1);
        
        if size(Polygon_intersection,2)~=0 % if the intersection is not empty...
            
            % Map of the neighnour
            BBox_ne = femregion_h.BBox(il,:);
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
            
            %Btilde11 = BJ_inv_ne(1,1)*BJ(1,1); Btilde22 = BJ_inv_ne(2,2)*BJ(2,2);
            %ttilde1 = BJ_inv_ne(1,1)*trans(1) + trans_inv_ne(1);
            %ttilde2 = BJ_inv_ne(2,2)*trans(2) + trans_inv_ne(2);
            Btilde11 = BJ_inv(1,1)*BJ_ne(1,1); Btilde22 = BJ_inv(2,2)*BJ_ne(2,2);
            ttilde1 = BJ_inv(1,1)*trans_ne(1) + trans_inv(1);
            ttilde2 = BJ_inv(2,2)*trans_ne(2) + trans_inv(2);
            
            nedges = length(Polygon_intersection.x);
            
            v = [Polygon_intersection.x(end:-1:1), Polygon_intersection.y(end:-1:1)]';
            
            % Define \Omega^ = Finv(\Omega)
            vhat = [];
            for i = 1:size(v,2)
                %vhat = [vhat, BJ_inv*v(:,i) + trans_inv];
                vhat = [vhat, BJ_inv_ne*v(:,i) + trans_inv_ne];
            end
            
            %Homoenenous mass matrix on \Omega^
            HM = [];
            for k = 0:2*p
                for l = 0:2*p-k
                    HM(k+1,l+1) = IntegralOverPolygon(1, k, l, vhat, toll);
                end
            end
            
            % Local coeff for neighbour shape functions
            for k = 0:p
                for j = k:p
                    
                    for l = k:j
                        %Ctilde1(j+1,k+1) = Ctilde1(j+1,k+1) + ...
                        %    P_coeff(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde11^k * ttilde1^(l-k);
                        %
                        %Ctilde2(j+1,k+1) = Ctilde2(j+1,k+1) + ...
                        %    P_coeff(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde22^k * ttilde2^(l-k);
                        
                        Ctilde1(j+1,k+1) = Ctilde1(j+1,k+1) + ...
                            P_coeff(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde11^k * ttilde1^(l-k);
                        
                        Ctilde2(j+1,k+1) = Ctilde2(j+1,k+1) + ...
                            P_coeff(j+1,l+1) * bin_coeff(l+1,k+1) * Btilde22^k * ttilde2^(l-k);
                    end
                    
                end
            end
            
            
            for idx = 0:p
                for idy = 0:p
                    
                    for k = 0:idx+idy
                        Cxy1 = 0; Cxy2 = 0;
                        for i = 0:min(k,idx)
                            if 0 <= k-i && k-i <= idy
                                %Cxy1 = Cxy1 + P_coeff(idx+1,i+1)*Ctilde1(idy+1,k-i+1);
                                %Cxy2 = Cxy2 + P_coeff(idx+1,i+1)*Ctilde2(idy+1,k-i+1);
                                
                                Cxy1 = Cxy1 + Ctilde1(idx+1,i+1)*P_coeff(idy+1,k-i+1);
                                Cxy2 = Cxy2 + Ctilde2(idx+1,i+1)*P_coeff(idy+1,k-i+1);
                            end
                        end
                        Ctilde1_x1_x2(idx+1,idy+1,k+1) = Cxy1;
                        Ctilde2_x1_x2(idx+1,idy+1,k+1) = Cxy2;
                    end
                    
                end
            end
            
            for idI = 1:femregion_h.nln
                
                % Phi_I(x,y) = phi_{ix}(x) * phi_{iy}(y)
                ix = local_enumeration_basis_functions(idI,1);
                iy = local_enumeration_basis_functions(idI,2);
                
                for idJ = 1:femregion_H.nln % and then consider the symmetry of M and V
                    
                    jx = local_enumeration_basis_functions(idJ,1);
                    jy = local_enumeration_basis_functions(idJ,2);
                    
                    Int = 0;
                    for k = 0:ix+jx
                        for l = 0:iy+jy
                            
                            %Int = Int + Ctilde1_x1_x2(ix+1,jx+1,k+1) * Ctilde2_x1_x2(iy+1,jy+1,l+1) ...
                            %    * Jdet * HM(k+1,l+1);
                            
                            Int = Int + Ctilde1_x1_x2(ix+1,jx+1,k+1) * Ctilde2_x1_x2(iy+1,jy+1,l+1) ...
                                * Jdet_ne * HM(k+1,l+1);
                            
                        end
                    end
                    
                    M_hH(index_h(idJ),index_H(idI)) = Int;
                   
                end
                
                
            end
            
            
            
        end
    end
    
end



end %end function