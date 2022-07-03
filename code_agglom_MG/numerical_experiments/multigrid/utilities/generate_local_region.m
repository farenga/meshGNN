function [region] = generate_local_region(Node,Element)

ne = length(Element);
elem_area = zeros(ne,1);
nedge = zeros(ne,1);
BBox = zeros(ne,4);

coords_element = cell(1,ne);
max_kb = cell(1,ne);

coord = Node;

for i = 1:length(Element)
    nedge(i) = length(Element{i});
    coords_element{i} = coord(Element{i},:);
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



region=struct('nedges', nedge',...
    'BBox',BBox,...
    'ne',length(Element),...
    'coord',coord);
region.coords_element = coords_element;
region.connectivity = Element;
region.area = elem_area;
region.max_kb = max_kb;


% Evaluate diameteer
for i = 1:region.ne
    coordi = region.coords_element{i};
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
    for j = 1:length(Element{i})
        if j<length(Element{i})
            edges{i} = [edges{i}; j j+1];
            edges_phys{i} = [edges_phys{i}; Element{i}(j) Element{i}(j+1)];
        else
            edges{i} = [edges{i}; j 1];
            edges_phys{i} = [edges_phys{i}; Element{i}(j) Element{i}(1)];
        end
    end
end
%region.h = max(H);
region.h = H;
region.edges = edges;
region.edges_phys = edges_phys;

