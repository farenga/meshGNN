function [u] = TWO_LEVELS_AS(A,f,z0,m1,m2,R_h,local_solvers,coarse_solver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MULTIGRID TWO LEVEL ALGORITHM WITH ASM  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(A)>2)
    X = 'Error: number of levels is different from 2.';
    disp(X);
    exit;
end

%Iniatialize
u = smoother_pcg(A{1}, z0, f, m1, local_solvers, coarse_solver);

%Coarse grid correction
%y = action_preconditioner(f - A{1}*u, local_solvers, coarse_solver, 'ADDITIVE', A{1});
y = f - A{1}*u;
r = R_h{1}'*(y);
e = A{2}\r;
u = u + R_h{1}*e;

%Post-smoothing
u = smoother_pcg(A{1}, u, f, m2, local_solvers, coarse_solver);


end % end function


