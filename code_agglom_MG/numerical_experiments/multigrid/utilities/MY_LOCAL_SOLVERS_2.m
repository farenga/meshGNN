%======================================================================
% Here we build the local solvers
%======================================================================

function [local_solvers] = MY_LOCAL_SOLVERS_2(subdomains_region, femregion, A)
% "A" and "femregion" are related to the fine space, which is usually the one
% with the subscript _h;
% "subdomains_region" is the one employed for local solvers, basically its
% elements are obtained by agglomeration of fine elements.

%disp('-------- BEGIN COMPUTING LOCAL SOLVER --------');

%===========================================================
% OUTPUT
local_solvers = struct('NSUB_i',{},...
    'Ai',{},...
    'Ri_T',{},...
    'index_dof',{});
%===========================================================

Nsub = subdomains_region.ne;

ndof = femregion.ndof;
Ne = femregion.ne; fine_element_loop = [1:Ne];
%ndof_i = femregion.nln;

for isub = 1:Nsub
    
    connect_subdom = subdomains_region.connectivity{isub};
    Xsub = subdomains_region.coord(connect_subdom,1);
    Ysub = subdomains_region.coord(connect_subdom,2);
    
    index = []; element_in_this_subdomain = [];
    for ie = fine_element_loop
        connect_element = femregion.connectivity{ie};
        xcoord = femregion.coord(connect_element,1); 
        ycoord = femregion.coord(connect_element,2);
        
        xm = sum(xcoord)/length(xcoord); ym = sum(ycoord)/length(ycoord);
        
        if inpolygon(xm,ym,Xsub,Ysub)
            element_in_this_subdomain = [element_in_this_subdomain, ie];
            index = [index; (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]'];
        end
    end
    %fine_element_loop = setdiff(fine_element_loop,element_in_this_subdomain);
    
    ndof_i = femregion.nln*length(element_in_this_subdomain);
    
    %==============================================
    % BUILD Ri_T
    %==============================================
    
    Ri_T = sparse(ndof,ndof_i);
    Ri_T(index,:) = speye(ndof_i);
    
    A_i =  transpose(Ri_T) * A * Ri_T;
    
    %==============================================
    % BUILD OUTPUT for the domain Omega_i
    %==============================================
    local_solvers(isub).NSUB_i = isub;
    local_solvers(isub).Ai = A_i;
    local_solvers(isub).index_dof = index;
    local_solvers(isub).Ri_T = Ri_T;
end


