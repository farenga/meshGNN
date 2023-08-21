function [x] = smoother_pcg(A, x0, b, max_it, local_solvers, coarse_solver)

x=x0;

r=b-A*x;

y = action_preconditioner(r, local_solvers, coarse_solver, 'ADDITIVE', A);
alpha=0;

for iter = 1:max_it                       % begin iteration
    rho = (r'*y);
    
    if ( iter > 1 )                       % direction vector
        beta = rho / rho_1;
        p = y + beta*p;
    else
        p = y;
    end
    
    q=A*p;
    alpha = rho / (p'*q );
    x = x + alpha * p;                    % update approximation vector
    r = r - alpha*q;                      % compute residual
    y = action_preconditioner(r, local_solvers, coarse_solver, 'ADDITIVE', A);
    rho_1 = rho;
end

end
