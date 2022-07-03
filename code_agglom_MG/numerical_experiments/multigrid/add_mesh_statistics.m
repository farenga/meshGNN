function region = add_mesh_statistics(region)
ne = region.ne;
Element = region.connectivity;
coord = region.coord;


elem_area = zeros(ne,1);
nedges = zeros(1,ne);
BBox = zeros(ne,4);
coords_element = cell(1,ne);
max_kb = cell(1,ne);

for i = 1:ne
    nedges(i) = length(Element{i});
    coords_element{i} = coord(Element{i},:);
    if ispolycw(coords_element{i}(:,1), coords_element{i}(:,2))
        [coords_element{i}(:,1), coords_element{i}(:,2)] = poly2ccw(coords_element{i}(:,1), coords_element{i}(:,2));
        region.connectivity{i} = fliplr(Element{i});
    end
    x_min = min( coords_element{i}(:,1));x_max = max( coords_element{i}(:,1));
    y_min = min( coords_element{i}(:,2));y_max = max( coords_element{i}(:,2));
    BBox(i,:)=[x_min x_max y_min y_max];
    
    P = polyshape(coords_element{i}(:,1),coords_element{i}(:,2),'KeepCollinearPoints',true);
    elem_area(i) = area(P);
    
    max_kb{i} = zeros(nedges(i),1);
    for j = 1:nedges(i)
        if j<nedges(i)
            j_plus = j+1;
        else
            j_plus = 1;
        end
        v1 = coords_element{i}(j,:);
        v2 = coords_element{i}(j_plus,:);
        
        for k = 1:nedges(i)
            if k~=j && k~=j_plus
                v3 = coords_element{i}(k,:);
                P_tria = polyshape([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)],'Simplify',false);
%                 P_int = intersect(P,P_tria);
%                 A_int = area(P_int);
                A_tria = area(P_tria);
                if isequal(subtract(P_tria,P),polyshape())
                    if A_tria>max_kb{i}(j)
                        max_kb{i}(j) = A_tria;
                    end
                end
               
            end
        end

        
    end
end



region.nvert = size(region.coord,1);
region.nedges = nedges;
region.BBox = BBox;
region.coords_element = coords_element;
region.area = elem_area;
region.max_kb = max_kb;
region.markers = 0;