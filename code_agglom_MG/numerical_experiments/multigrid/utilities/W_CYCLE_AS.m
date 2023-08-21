function [u] = W_CYCLE_AS(total_levels,lev,A,f,z0,m1,m2,R_h,local_solvers,coarse_solver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   MULTIGRID TWO LEVEL ALGORITHM   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if(length(A) > total_levels)
%     error('Error: number of levels is different from the number of total levels.');
% end

%Iniatialize
u = z0;

if(lev == total_levels)

    u = A{lev}\f;

else
    
    
    % Pre-smoothing
    u = smoother_pcg(A{lev}, z0, f, m1, local_solvers{lev}, coarse_solver{lev});
    
    
    % Coarse grid correction
    %y = action_preconditioner(f - A{lev}*u, local_solvers{lev}, coarse_solver{lev}, 'ADDITIVE', A{lev});
    y = f - A{lev}*u;
    %y = A{lev}*u - f;
    r = (R_h{lev}')*y;
    start = zeros(length(r),1);
    e = W_CYCLE_AS(total_levels,lev+1,A,r,start,m1,m2,R_h,local_solvers,coarse_solver);
    e = W_CYCLE_AS(total_levels,lev+1,A,r,e,m1,m2,R_h,local_solvers,coarse_solver);
    u = u + R_h{lev}*e;
    
    % Post-smoothing
    u = smoother_pcg(A{lev}, u, f, m2, local_solvers{lev}, coarse_solver{lev});

    
end