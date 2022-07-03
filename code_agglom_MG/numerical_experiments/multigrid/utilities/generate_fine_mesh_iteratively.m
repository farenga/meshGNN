function [Tau_fine] = generate_fine_mesh_iteratively(Tau_coarse,Nit)

Dati = dati_2;
Dati = dati_2;
Dati.type_mesh = 'QUAD';
Dati.method = 'LDG';
Dati.basis = 'modal';
[Tau_cartesian] = generate_mesh_2(Dati,3);

Tau_cartesian.coords(1,:) = 2*Tau_cartesian(1,:)-1;
for ie=1:Tau_cartesian.ne
    Tau_cartesian{ie}(1,:) = 2*Tau_cartesian{ie}(1,:)-1;
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
for ie = 1:Tau_coarse.ne
    
    P1.x    = [Tau_coarse.coords_element{ie}(:,1);Tau_coarse.coords_element{ie}(1,1)];
    P1.y    = [Tau_coarse.coords_element{ie}(:,2);Tau_coarse.coords_element{ie}(1,2)];
    P1.hole = 0;
    
    xH_min = Tau_coarse.BBox(ie,1); xH_max = Tau_coarse.BBox(ie,2);
    yH_min = Tau_coarse.BBox(ie,3); yH_max = Tau_coarse.BBox(ie,4);
    
    xH_half = (xH_min + xH_max)/2; yH_half = (yH_min + yH_max)/2;
    
    
    % Element 1: intersection with lower right
    P2.x = [xH_half xH_max xH_max xH_half xH_half];
    P2.y = [yH_min yH_min yH_half yH_half yH_min];
    P1.hole = 0;
    
    P = PolygonClip(P1,P2,1);
    if size(P,2)~=0
        ne = ne + 1;
        nedges(ne) = length(P.x);
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
    
    
    
    
    
    
    % Element 2: intersection with upper right
    P2.x = [xH_half xH_max xH_max xH_half xH_half];
    P2.y = [yH_half yH_half yH_max yH_max yH_half];
    P1.hole = 0;
    
    P = PolygonClip(P1,P2,1);
    if size(P,2)~=0
        ne = ne + 1;
        nedges(ne) = length(P.x);
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
    
    
    
    
    
    % Element 3: intersection with upper left
    P2.x = [xH_min xH_half xH_half xH_min xH_min];
    P2.y = [yH_half yH_half yH_max yH_max yH_half];
    P1.hole = 0;
    
    P = PolygonClip(P1,P2,1);
    if size(P,2)~=0
        ne = ne + 1;
        nedges(ne) = length(P.x);
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
    
    
    
    
    
    
    % Element 4: intersection with lower left
    P2.x = [xH_min xH_half xH_half xH_min xH_min];
    P2.y = [yH_min yH_min yH_half yH_half yH_min];
    P1.hole = 0;
    
    P = PolygonClip(P1,P2,1);
    if size(P,2)~=0
        ne = ne + 1;
        nedges(ne) = length(P.x);
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


Tau_fine.nedges = nedges;
Tau_fine.BBox = BBox;
Tau_fine.ne = ne;
Tau_fine.coord = coord;
Tau_fine.coords_element = coords_element;
Tau_fine.connectivity = connectivity;

