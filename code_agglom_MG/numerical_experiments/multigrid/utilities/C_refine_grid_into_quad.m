function [refined_region]=C_refine_grid_into_quad(region, marked_elem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine marked elements of a grid
% Input: region: structure with pointers
%                coord->  coordinated of vertices
%                ne   ->  number of elements
%                connectivity ->  connectivity vectors for each mesh
%                nvert   ->  number of vertexes
%        marked_elem: vector of size region.ne. It contains 1 if the
%        element has to be refined and 0 otherwise
% Output: refined_region (structure with the same pointers as above)
% Author: Paola F. Antonietti

%
% Giorgio Pennesi modified something 22/03/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cd ..
coord=region.coord;
 
if nargin ==1
    marked_elem=ones(region.ne,1);
end
 
% add centers and midpoints of edges
% loop over the elements
tk=1;
for ie=1:region.ne
    
    if marked_elem(ie)==1 %refine
        % get the connectivity of the element
        %[index]=C_get_index_element(region.connectivity, ie);
        index = region.connectivity{ie};
        
        % midpoints
        [mp]=C_get_midpoints(region.coord(index,:));
                
        % centroid
        BB=sum(region.coord(index,:))./(length(region.coord(index,1)));
        
        % check for duplications of midpoints
        newcoord=[];
        nindex=[];
        t=1;
        for k=1:size(mp, 1)
            tmp=find(mp(k,1)==coord(:,1) &  mp(k,2)==coord(:,2));
            if isempty(tmp)
                newcoord=[newcoord; mp(k,:)];
                nindex(k)=[(size(coord,1))+t];
                t=t+1;
            else
                nindex(k)=tmp;
            end
        end
        nindex=[nindex,size(coord,1)+size(newcoord,1)+1];
        
        % add coordinates
        coord=[coord;newcoord;BB];
        
        % 0 for midpoints,1 for center of mass 
        for k=1:length(index)
            v1=index(k);
            v2=nindex(k);
            v3=nindex(end);
            if (k==1)
                v4=nindex(end-1);
            else
                v4=nindex(k-1);
            end
            connectivity{tk}=[v1;v2;v3;v4];
            tag(tk) = region.tag(ie);
            tk=tk+1;
        end
    end
end
 
new_dof=coord(size(region.coord,1)+1:end,:);
% update connectivity unmarked elements
if prod(marked_elem)==0
    unmarked_elemts=find(marked_elem==0);
    for k=1:length(unmarked_elemts)
        [index]=C_get_index_element(region.connectivity, unmarked_elemts(k));     
        v=[index(1)];
        %loop over vertexes
        for j=1:length(index)
            iin=index(j);
            vin=coord(iin,:);
            if j==length(index)
                ifin=index(1);
                vfin=coord(index(1),:);
            else
                ifin=index(j+1);
                vfin=coord(index(j+1),:);
            end
            vmid=(vin+vfin)./2;
            flag=find(vmid(1)==new_dof(:,1) & vmid(2)==new_dof(:,2));
            if isempty(flag)
                v=[v,ifin];
            else
                v=[v,flag+region.nvert,ifin];
            end
        end
       connectivity{tk}=v(1:end-1)';
       tk=tk+1;
    end
end
    
 
% update region
refined_region.ne=size(connectivity,2);
refined_region.nvert=size(coord,1);
refined_region.coord=coord;
refined_region.connectivity=connectivity;


ne = refined_region.ne;
elem_area = zeros(ne,1);
nedge = zeros(ne,1);
BBox = zeros(ne,4);

coords_element = cell(1,ne);
max_kb = cell(1,ne);

%coord = Node;

for i = 1:ne
    nedge(i) = length(refined_region.connectivity{i});
    coords_element{i} = coord(refined_region.connectivity{i},:);
    x_min = min( coords_element{i}(:,1));x_max = max( coords_element{i}(:,1));
    y_min = min( coords_element{i}(:,2));y_max = max( coords_element{i}(:,2));
    BBox(i,:)=[x_min x_max y_min y_max];
    
    elem_area(i) = polyarea(coords_element{i}(:,1),coords_element{i}(:,2));
    max_kb{i} = zeros(nedge(i),1);
    for j = 1:nedge(i)
        if j<nedge(i)
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(j+1,:);
        else
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(1,:);
        end
        for k = 1:nedge(i)
            if k~=j && k~=(j+1)
                v3 = coords_element{i}(k,:);
                [x_tria,y_tria]=poly2cw([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)]);
                [x1,y1] = polybool('intersection',coords_element{i}(end:-1:1,1),coords_element{i}(end:-1:1,2),x_tria,y_tria);
                area = polyarea(x_tria,y_tria);
                if 1-any(isnan(x1)) && abs(polyarea(x1,y1)- area)<1e-13
                    if area>max_kb{i}(j)
                        max_kb{i}(j) = area;
                    end
                end
            end
        end
        
    end
    
end



refined_region.nedges= nedge;
refined_region.BBox=BBox;
    
refined_region.coords_element = coords_element;
%refined_region.connectivity = Element;
refined_region.area = elem_area;
refined_region.max_kb = max_kb;


% Evaluate diameteer
for i = 1:refined_region.ne
    coordi = refined_region.coords_element{i};
    cc = 1;
    for j = 1:size(coordi,1)-1
        for k = j+1:size(coordi,1)
            distance(cc) = sqrt((coordi(j,1)-...
                coordi(k,1)).^2+(coordi(j,2)-coordi(k,2)).^2);
            cc = cc+1;
        end
    end
    H(i) = max(distance);
    
    edges{i}=[];
    edges_phys{i}=[];
    for j = 1:length(refined_region.connectivity{i})
        if j<length(refined_region.connectivity{i})
            edges{i} = [edges{i}; j j+1];
            edges_phys{i} = [edges_phys{i}; refined_region.connectivity{i}(j) refined_region.connectivity{i}(j+1)];
        else
            edges{i} = [edges{i}; j 1];
            edges_phys{i} = [edges_phys{i}; refined_region.connectivity{i}(j) refined_region.connectivity{i}(1)];
        end
    end
end
%refined_region.hmax = max(H);
refined_region.h = H;
refined_region.edges = edges;
refined_region.edges_phys = edges_phys;

refined_region.tag = tag;








    
 
% OPTIONAL
% % get baricenters
% for ie=1:refined_region.ne
%     [index]=C_get_index_element(refined_region.connectivity, ie);
%     BB=sum(refined_region.coord(index,:))./(length(refined_region.coord(index,1)));
%     refined_region.baricenters(ie,:)=BB;
% end

