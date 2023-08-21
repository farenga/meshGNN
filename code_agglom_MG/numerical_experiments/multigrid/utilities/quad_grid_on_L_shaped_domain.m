function [Tau_fine] = quad_grid_on_L_shaped_domain(N)

DatiHere = dati_2;
DatiHere.type_mesh = 'QUAD';
DatiHere.method = 'LDG';
DatiHere.basis = 'modal';
[Tau_cartesian] = generate_mesh_2(DatiHere,N);

Tau_cartesian.BBox(:,1) = 2*Tau_cartesian.BBox(:,1);
Tau_cartesian.BBox(:,2) = 2*Tau_cartesian.BBox(:,2);
Tau_cartesian.BBox(:,3) = 2*Tau_cartesian.BBox(:,3);
Tau_cartesian.BBox(:,4) = 2*Tau_cartesian.BBox(:,4);

Tau_cartesian.coord(:,1) = 2*Tau_cartesian.coord(:,1);
Tau_cartesian.coord(:,2) = 2*Tau_cartesian.coord(:,2);

for ie=1:Tau_cartesian.ne
    Tau_cartesian.coords_element{ie}(:,1) = 2*Tau_cartesian.coords_element{ie}(:,1);
    Tau_cartesian.coords_element{ie}(:,2) = 2*Tau_cartesian.coords_element{ie}(:,2);
end



nedges = [];
BBox = [];
ne = 0;
coord = [];
coords_element = {};
connectivity = {};
area = {};
max_kb = {};
h = [];
for ie = 1:Tau_cartesian.ne
    
    P1.x    = [Tau_cartesian.coords_element{ie}(:,1);Tau_cartesian.coords_element{ie}(1,1)];
    P1.y    = [Tau_cartesian.coords_element{ie}(:,2);Tau_cartesian.coords_element{ie}(1,2)];
    P1.hole = 0;
    
    L_shaped.x    = [0 2 2 1 1 0 0]';
    L_shaped.y    = [0 0 1 1 2 2 0]';
    L_shaped.hole = 0;
    
    Pintersect = PolygonClip(P1,L_shaped,1);
    if(size(Pintersect,2)~=0)
        for i=1:size(Pintersect,2)
            P = Pintersect(i);
            if size(P,2)~=0
                P.x = P.x(end:-1:1);
                P.y = P.y(end:-1:1);
                ne = ne + 1;
                nedges(ne) = size(P.x,1);
                BBox(ne,1) = min(P.x); BBox(ne,2) = max(P.x);
                BBox(ne,3) = min(P.y); BBox(ne,4) = max(P.y);
                coords_element{ne} = [P.x P.y];
                
                connect_element = [];
                for i=1:nedges(ne)
                    xi = P.x(i); yi = P.y(i);
                    
                    j=1;
                    while j<=size(coord,1)
                        if ((xi-coord(j,1) == 0) && (yi-coord(j,2)==0))
                            connect_element = [connect_element j];
                            j = size(coord,1) + 10;
                        else
                            j = j+1;
                        end
                    end
                    j = j-1;
                    if j==size(coord,1)
                        coord = [coord' [xi yi]']';
                        connect_element = [connect_element size(coord,1)];
                    end
                    
                end
                connectivity{ne} = connect_element;
                
            end
        end
        
    end
    
    
end


Tau_fine.nedges = nedges;
Tau_fine.BBox = BBox;
Tau_fine.ne = ne;
Tau_fine.coord = coord;
Tau_fine.coords_element = coords_element;
Tau_fine.connectivity = connectivity;



% Evaluate diameteer
for i = 1:Tau_fine.ne
    coordi = Tau_fine.coords_element{i};
    cc = 1;
    for j = 1:size(coordi,1)-1
        for k = j+1:size(coordi,1)
            distance(cc) = sqrt((coordi(j,1)-...
                coordi(k,1)).^2+(coordi(j,2)-coordi(k,2)).^2);
            cc = cc+1;
        end
    end
    H(i) = max(distance);
end
Tau_fine.H = H;




max_kb = cell(1,ne);
for i = 1:ne
    elem_area(i) = polyarea(coords_element{i}(:,1),coords_element{i}(:,2));
    max_kb{i} = zeros(nedges(i),1);
    for j = 1:nedges(i)
        if j<nedges(i)
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(j+1,:);
        else
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(1,:);
        end
        for k = 1:nedges(i)
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

Tau_fine.area = area;
Tau_fine.max_kb = max_kb;

