%======================================================================
% Here we build the local solvers
%======================================================================

function [local_solvers] = MY_LOCAL_SOLVERS(femregion, A)

%disp('-------- BEGIN COMPUTING LOCAL SOLVER --------');

%===========================================================
% OUTPUT
local_solvers = struct('NSUB_i',{},...
                   'Ai',{},...
                   'Ri_T',{},...
                   'index_dof',{});
%===========================================================

ndof = femregion.ndof;
Ne = femregion.ne;
ndof_i = femregion.nln;

for ie = 1:Ne
    
    index = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    
    %==============================================
    % BUILD Ri_T
    %==============================================
    
    Ri_T = sparse(ndof,ndof_i);
    Ri_T(index,:) = speye(ndof_i);
    
    A_i =  transpose(Ri_T) * A * Ri_T;
    
    %==============================================
    % BUILD OUTPUT for the domain Omega_i
    %==============================================
    local_solvers(ie).NSUB_i = ie;
    local_solvers(ie).Ai = A_i;
    local_solvers(ie).index_dof = index;
    local_solvers(ie).Ri_T = Ri_T;    
end


