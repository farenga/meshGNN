%--------------------------------------------------------------------
% PURPOSE:
%
% This routine contains all the information  for running main, i.e.,
%
% DATI.method        -> (string) 'SIPG'
% DATI.domain        -> (2x2 matrix, real)  vertices of the domain 
%                       (supposed to be a rectangle)            
% DATI.exact_sol     -> (string) analytical expression of the exact solution
% DATI.source        -> (string) analytical expression of source term
% DATI.grad_exact_1  -> (string) analytical expression of first derivative 
%                       (with respect to the x-variable) of the exact solution  
% DATI.grad_exact_2  -> (string) analytical expression of first derivative 
%                       (with respect to the y-variable) of the exact solution 
% DATI.fem           -> (string) finite element space ['P1', 'P2', 'P3', 'P4']
% DATI.penalty_coeff -> (real) penalty pearameter
% DATI.nqn           -> (integer) number of 1D Gauss-Ledendre quadrature nodes in 1 
%                       dimension [1,2,3,4,5,6,.....]
% DATI.type_mesh     ->  (string) type of mesh to be used [strctured triangular mesh 'TS' 
%                        or unstructured triangular mesh 'TU']
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [DATI]=dati

DATI= struct(  'domain',           [0,1;0,1],...              
               'exact_sol',        '(x.^2-x).*(y.^2-y).*exp(x)',...    
               'source',           '-(x.^2+3*x).*exp(x).*(y.^2-y)-2*(x.^2-x).*exp(x)',...  
               'grad_exact_1',     '(x.^2+x-1).*(y.^2-y).*exp(x)',...    
               'grad_exact_2',     '(x.^2-x).*(2*y-1).*exp(x)',...    
               'fem',               1,...   
               'penalty_coeff',     10,... 
               'nqn',               3);

           
           