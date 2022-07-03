function [x, error, iter, flag, condi, esti] = my_pcg(A, x0, b, max_it, tol,local_solvers, coarse_solver)

%
%
%
% Conjugate Gradient method with preconditioning.
% no reorthogonalization.
%
% input   A        REAL la matrice in questione...
%         x0       REAL il initial guess(non so come dirlo in italiano!)
%         b        REAL il lato destro della equazione
%         max_it   INTEGER numero massimo d'iterazioni
%         tol      REAL la tolerance
%
% output  x        REAL la soluzione
%         error    REAL il errore in norma L2
%         iter     INTEGER le iterazione fatte
%         flag     INTEGER: 0 ... tutto aposto OK
%                           1 ... Non CONVERGE
%         condi       REAL estimated condition number
%         esti      real.. estima degli autovalori..
%
% OCCHIO: devi creare una funzione aparte che si chiame preac
%tipo...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y=preac(x,rpar)
%in preac gia hai messo tutto quello di cui hai bisogno..
%e poi ..
%
% y=R\(R'\x);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%

flag = 0;                                 % initialization
iter = 0;
condi = 1;
esti(1) = 1;

x = x0;

r = b-A*x;
bnrm2 = norm(b);

if  ( bnrm2 == 0.0 ) 
    bnrm2 = 1.0; 
end

%y = action_preconditioner(b, local_solvers, coarse_solver);
y = action_preconditioner(b, local_solvers, coarse_solver, 'ADDITIVE', A);
error = norm( y ) / bnrm2;
%fprintf(' PCG residual(%d) = %g\n',iter,error)
if ( error < tol ) return, end

alpha = 0;

delta = []; eta = [];

for iter = 1:max_it                       % begin iteration
    
    rho = (r'*y);
    
    if ( iter > 1 )                       % direction vector
        beta = rho / rho_1;
        p = y + beta*p;
    else
        p = y;
    end
    
    q=A*p;
    alphaold = alpha;
    alpha = rho / (p'*q );
    x = x + alpha * p;                    % update approximation vector
    
    r = r - alpha*q;                      % compute residual
    y = action_preconditioner(r, local_solvers, coarse_solver, 'ADDITIVE', A);
    
    error = norm(r)/bnrm2;            % comprobo la convergence
    
    
    % This is for condition number estimates
    if ( iter > 1 )                       % direction vector
        delta = [delta, (1/alpha) + (beta/alphaold)];
        eta = [eta, sqrt(beta)/alphaold];
    else
        delta = [delta, 1/alpha];
    end
    
    %X = ['    ',num2str(error)]; disp(X);
    if ( error <= tol )
        %fprintf(' PCG residual(%d) = %g\n',iter,error)
        break;
    end
    rho_1 = rho;
    
end

if ( error > tol ) 
    flag = 1; 
end  % no convergence

T = diag(delta) - diag(eta,-1) - diag(eta,+1);
esti = eig(T);
condi = max(esti)/min(esti);

end