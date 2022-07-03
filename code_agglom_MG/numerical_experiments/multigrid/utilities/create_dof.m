function [femregion]=create_dof(Dati,region)

ne=region.ne;
nedge=region.nedges;
degree = Dati.fem;

nln=0.5.*(degree+1).*(degree+2);

femregion=struct('fem',Dati.fem,...
    'nedges', nedge,...
    'nln',nln,...
    'ndof',nln*ne,...
    'ne',ne,...
    'nqn',Dati.nqn,...
    'coord',region.coord,...
    'BBox',region.BBox);
    
femregion.coords_element = region.coords_element;
femregion.connectivity=region.connectivity;
femregion.area = region.area;
femregion.max_kb = region.max_kb;


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
end
region.h = H;
femregion.h = region.h;

