function [midpoints]=C_get_midpoints(coord_ie)
 
coord_ie_tmp=[coord_ie;coord_ie(1,:)];
 
for k=1:length(coord_ie(:,1))
        midpoints(k,:)=[coord_ie_tmp(k,1) +  coord_ie_tmp(k+1,1); coord_ie_tmp(k,2) +  coord_ie_tmp(k+1,2)]./2;  % corresponding normal
end