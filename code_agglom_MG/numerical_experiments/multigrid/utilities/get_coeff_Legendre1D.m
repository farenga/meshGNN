function [P_coeff,dP_coeff] = get_coeff_Legendre1D(deg)
%
% [P,dP] = GetCoeffLegendre1DAndDerivatives(deg)
% 
% Generate two matrices P and dP of size (deg+1) \times (deg+1), containing
% the coefficient lists of the Legendre 1D polynomials up to degree deg, 
% and their derivatives, respectively:
%  - P(1,:) contains the list of the Legendre polynomial of degree 0;
%  - P(k,:)    "     "    "   "   "     "        "       "   "    k-1;
%  - dP(k,:)   "     "    "   "   "  derivative of the Legendre polynomial
%         of degree k-1;
%
% Symbolic evaluation of Legendre Polynomials is used in oreder to write
% each of them them as a sum of homogeneous polynomial functions. Instead 
% the recoursive definition, we use the Rodrigue's Formula:
%             1    d^n
% Pn(x) = -------- ---- (x^2 - 1)^n
%         2^n * n! dx^n
%
% and the normalization Pn(x) = Pn(x)/||Pn||_{L2}.
%
% The goal is to write 
%               Pn(x) = C0 + C1*x + C2*x^2 + ... + Cn*x^n,
% and hence collect the coefficient vector 
%               P(n+1,:) = [C0, C1, ... , Cn].
% The function gives also the coefficient list of the polynomial's 
% derivatives, that is 
%               dPn
%               --- (x)/dx = C1 + 2*C2*x + ... + n*Cn*x^(n-1) 
%                dx
% will be represented by the coefficient vector 
%               dP(n+1,:) = [C1, 2*C1, ... , n*Cn, 0].
%
% Remark: the 0 at the end of dP(n+1,:) is because we want P(n+1,;) and
% dP(n+1,;) to be the same length, because it is helpful for the
% implementation.

x = sym('x');
P_coeff = zeros(deg+1,deg+1);
dP_coeff = zeros(deg+1,deg+1);
for p = 0:deg
    E = [];
    C = 1/((2^p) * factorial(p));
    poly = (x^2 - 1)^p;
    for j=1:p                    % evaluation of   d^n
        poly = diff(poly,x);     %                 ---- (x^2 - 1)^n
    end                          %                 dx^n
    poly = C * poly;
    E = collect(poly);
    P{p+1} = sym2poly(E);
    P{p+1} = sqrt((2*p+1)/2)*P{p+1}(end:-1:1);
    
    P_coeff(p+1,[1:length(P{p+1})]) = P{p+1}; 
    
    
    % Derivatives
    if p == 0
        dP{1}(1) = 0; dP{1}(2) = 0;
    else
        for i = 1:p
            dP{p+1}(i) = P{p+1}(i+1)*i;
        end
        dP{p+1}(p+1) = 0;
        dP{p+1}(p+2) = 0;
    end
    dP_coeff(p+1,[1:length(dP{p+1})]) = dP{p+1};
end


end