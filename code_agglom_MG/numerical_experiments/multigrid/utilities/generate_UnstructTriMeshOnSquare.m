function [region]=generate_UnstructTriMeshOnSquare(Dati,hmax)

nedge=3;
[P,E,T]=initmesh(@squareg,'Hmax',hmax);
v = [Dati.domain(1,1) Dati.domain(2,1); ...
     Dati.domain(1,2) Dati.domain(2,1); ...
     Dati.domain(1,2) Dati.domain(2,2); ...
     Dati.domain(1,1) Dati.domain(2,2)];
x1 = min(v(:,1)); x2 = max(v(:,1));
y1 = min(v(:,2)); y2 = max(v(:,2));
trans = [x1+x2, y1+y2]/2;
Lx = abs(x2 - x1); Ly = abs(y2 - y1);
if Lx ~= Ly
    fprintf('Warning: the mesh is not uniform.\n');
end

P(1,:) = Lx*0.5*P(1,:) + trans(1);
P(2,:) = Ly*0.5*P(2,:) + trans(2);
%pdemesh(P,E,T);

T(end,:) = [];
coord = P;
connectivity = T;
clear P; clear T;

ne=size(connectivity,2);
h = sqrt(polyarea(v(:,1),v(:,2))/ne);

connectivity = connectivity(1:nedge,:);

coords_element = []; % coordinates of the elements (counted with their multiplicity)
BBox = zeros(ne,4);
elem_area = zeros(ne,1);
max_kb = cell(1,ne);

%coords_element = cell(1,ne);

%coord = Node;
nedges = zeros(ne,1);
for ie=1:ne
    nedges(ie) = nedge;
    %loc_coord = [];
    for k=1:nedge
        coords_element = [coords_element; coord(1,connectivity(k,ie)),coord(2,connectivity(k,ie))];
        %loc_coord = [loc_coord; coord(1,connectivity(k,ie)),coord(2,connectivity(k,ie))];
    end
    %coords_element = [coords_element; coord(1,connectivity(1,ie)),coord(2,connectivity(1,ie))];
    %coords_element{ie} = loc_coord; plot([coords_element(1:4,1);coords_element(1,1)],[coords_element(1:4,2);coords_element(1,2)]); hold on
    %coords_element{ie} = coord(:,connectivity(:,ie)); plot([coords_element(5:8,1);coords_element(5,1)],[coords_element(5:8,2);coords_element(5,2)]);
    x_min = min( coords_element((ie-1)*nedge+1:ie*nedge,1) ); x_max = max( coords_element((ie-1)*nedge+1:ie*nedge,1) );
    y_min = min( coords_element((ie-1)*nedge+1:ie*nedge,2) ); y_max = max( coords_element((ie-1)*nedge+1:ie*nedge,2) );
    %x_min = min( loc_coord(:,1)); x_max = max( loc_coord(:,1));
    %y_min = min( loc_coord(:,2)); y_max = max( loc_coord(:,2));
    BBox(ie,:)=[x_min x_max y_min y_max];
    
    elem_area(ie) = polyarea(coords_element((ie-1)*nedge+1:ie*nedge,1),coords_element((ie-1)*nedge+1:ie*nedge,2));
    %elem_area(ie) = polyarea(loc_coord(:,1),loc_coord(:,2));
    max_kb{ie} = zeros(nedge,1);
    
    for j = 1:nedge
        if j<nedge
            v1 = coords_element((ie-1)*nedge+j,:); v2 = coords_element((ie-1)*nedge+(j+1),:);
            %v1 = loc_coord(j,:); v2 = loc_coord(j+1,:);
        else
            v1 = coords_element((ie-1)*nedge+j,:); v2 = coords_element((ie-1)*nedge+1,:);
            %v1 = loc_coord(j,:); v2 = loc_coord(1,:);
        end
        for k = 1:nedge
            if k~=j && k~=(j+1)
                v3 = coords_element((ie-1)*nedge+k,:);
                %v3 = loc_coord(k,:);
                [x_tria,y_tria]=poly2cw([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)]);
                local_coords = coords_element((ie-1)*nedge+1:ie*nedge,:);
                [x1,y1] = polybool('intersection',local_coords(end:-1:1,1),local_coords(end:-1:1,2),x_tria,y_tria);
                %[x1,y1] = polybool('intersection',loc_coord(end:-1:1,1),loc_coord(end:-1:1,2),x_tria,y_tria);
                area = polyarea(x_tria,y_tria);
                if 1-any(isnan(x1)) && abs(polyarea(x1,y1)- area)<1e-13
                    if area>max_kb{ie}(j)
                        max_kb{ie}(j) = area;
                    end
                end
            end
        end
    end
end

Element = cell(1,ne); coords_element_cell = cell(1,ne);
for ie=1:ne
    Element{ie} = connectivity(:,ie);
    coords_element_cell{ie} = coords_element((ie-1)*nedge+1:ie*nedge,:);
end

region=struct('dim',2,...
    'domain',Dati.domain,...
    'nedges', nedges,...
    'h',h,...
    'nvert',size(coord,2),...
    'ne',ne,...
    'coord',coord',...
    'connectivity',connectivity,...
    'coords_element',coords_element);

%Added by GP
region.BBox = BBox;
region.area = elem_area;
region.max_kb = max_kb;
region.connectivity = Element;
region.coords_element = coords_element_cell;

end