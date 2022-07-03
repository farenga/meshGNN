function [x, error, iter, flag, condi, esti] = my_cg(A, x0, b, max_it, tol)

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
  condi=1;
  esti(1)=1;

  x=x0;

  r=b-A*x;
  bnrm2=norm(b);

  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
  y = r;
  error = norm( y ) / bnrm2;
  %fprintf(' PCG residual(%d) = %g\n',iter,error)
  if ( error < tol ) return, end


  alpha=0;
  
  
for iter = 1:max_it                       % begin iteration
     rho = (r'*y);
    
     if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = y + beta*p;
     else
        p = y;
     end

     q=A*p;
     alphaold=alpha;
     alpha = rho / (p'*q );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual
     y = r;
     
     error = norm(r)/bnrm2;            % comprobo la convergence
     %fprintf(' PCG residual(%d) = %g\n',iter,error)

     if ( iter > 1 )                       % matrice tridiagonal 
       tri(iter,iter)=(1/alpha)+(beta/alphaold);
       tri(iter,iter-1)=-sqrt(beta)/alphaold;
       tri(iter-1,iter)=tri(iter,iter-1);
     else
       tri(iter,iter)=(1/alpha);
     end 
     
     if ( error <= tol )
       esti=eig(tri);
       condi=max(esti)/min(esti);
       fprintf(' PCG residual(%d) = %g\n',iter,error)
       break;
   end 

     rho_1 = rho;
  end
  esti=eig(tri);
  condi=max(esti)/min(esti);
  if ( error > tol ), flag = 1; end  % no convergence

