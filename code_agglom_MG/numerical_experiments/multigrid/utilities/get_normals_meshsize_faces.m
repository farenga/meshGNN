%--------------------------------------------------------------------
% PURPOSE:
% This routine computes, for each face of a given element, the diameter 
% of the face and the corresponding unit outward normal
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [normal,meshsize]=get_normals_meshsize_faces(loc_coord)

nfaces=length(loc_coord);
normal = zeros(2, nfaces);
meshsize= zeros(nfaces,1);

loc_coord_tmp = [loc_coord; loc_coord(1,1), loc_coord(1,2)];

for i = 1:nfaces              
    nn=[-loc_coord_tmp(i,2) +  loc_coord_tmp(i+1,2); loc_coord_tmp(i,1) -  loc_coord_tmp(i+1,1)];  % corresponding normal
    nn=nn/norm(nn);                            % unit normal 
    normal(:,i)=nn;
    meshsize(i,1) = norm([loc_coord_tmp(i+1,1)- loc_coord_tmp(i,1) ; loc_coord_tmp(i+1,2)- loc_coord_tmp(i,2)]);
    clear nn
end

%GP: if F is a face then all [x,y] \in F satisfy <A,[x,y]> = b, 
% where A = [a1,a2] and b is a number.
% if the element has a counterclockwise node list, and F is the line
% segment between [x1,y1] and [x2,y2], then it is possoble to choose
%                 1                 | y2 - y1 |
% A = ----------------------------- |         |   
%     sqrt((y2-y1)**2 + (x1-x2)**2) | x1 - x2 |
%
% b = [(y2-y1)*x1 + (x1-x2)*y1]/sqrt((y2-y1)**2 + (x1-x2)**2) 

