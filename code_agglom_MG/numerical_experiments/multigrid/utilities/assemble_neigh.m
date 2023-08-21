%--------------------------------------------------------------------
% PURPOSE:
%
% This routine assembles the local matrices corresponding
% to the neighbours of a given element.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [M]= assemble_neigh(M,row,neight,M1,nln,n_edge)

for iedg=1:n_edge
    if neight(iedg)~= -1
        j=(neight(iedg)-1)*nln*ones(nln,1) + [1:nln]';
        M(row,j)=M(row,j)+M1(:,:,iedg);
    end
end
