%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Modification of PolyMesher      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [region] = generate_mesh_for_L_shaped(Domain,NElem,MaxIter,P)
if ~exist('P','var'), P=PolyMshr_RndPtSet(NElem,Domain); end
NElem = size(P,1);
Tol=5e-3; It=0; Err=1; c=1.5;%c=1.5;
BdBox = Domain('BdBox');
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Pc = P;
while(It<=MaxIter && Err>Tol)
    Alpha = c*sqrt(Area/NElem);
    P = Pc;                                       %Lloyd's update
    R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha);   %Generate the reflections
    [Node,Element] = voronoin([P;R_P],{'Qbb','Qz'});           %Construct Voronoi diagram
    [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem);
    Area = sum(abs(A));
    Err = sqrt(sum((A.^2).*sum((Pc-P).*(Pc-P),2)))*NElem/Area^1.5;
    %   fprintf('It: %3d   Error: %1.3e\n',It,Err); It=It+1;
end
[Node,Element] = PolyMshr_ExtrNds(NElem,Node,Element);  %Extract node list
[Node,Element] = PolyMshr_CllpsEdgs(Node,Element,0.1);  %Remove small edges
[Node,Element] = PolyMshr_RsqsNds(Node,Element);        %Reoder Nodes
%BC=Domain('BC',Node); Supp=BC{1}; Load=BC{2};           %Recover BC arrays

%PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load);

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


%per capire quali sono le coordinate degli elementi dell'elemento
%coords_element{length(Element)}

%%% POLYMESHER FUNCTIONS %%%%

%------------------------------------------------- GENERATE RANDOM POINTSET
% GP modified this, 17/03/2019
function P = PolyMshr_RndPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;

%fixed seeds
% P(1,1) = 0.75; P(1,2) = 1.25; Ctr=Ctr+1;
% P(2,1) = 0.75; P(2,2) = 0.75; Ctr=Ctr+1;
% P(3,1) = 1.25; P(3,2) = 0.75; Ctr=Ctr+1;
help = NElem;
NElem = NElem-3;
while Ctr<NElem
    Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
    Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
    I = find( sqrt( (Y(:,1)-1).^2 + (Y(:,2)-1).^2 ) > 0.3 );
    Y = Y(I,:);
    d = Domain('Dist',Y);
    I = find(d(:,end)<0);                 %Index of seeds inside the domain
    NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
    P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
    Ctr = Ctr+NumAdded;
end

% [val,index1] = min(sqrt( (P(:,1)-1).^2 + (P(:,2)-1).^2 ));
% x1 = P(index1,1); y1 = P(index1,2);
% P(index1,1) = 0; P(index1,2) = 0;
%
% [val,index2] = min(sqrt( (P(:,1)-1).^2 + (P(:,2)-1).^2 ));
% x2 = P(index2,1); y2 = P(index2,2);
% P(index2,1) = 0; P(index2,2) = 0;
%
% [val,index3] = min(sqrt( (P(:,1)-1).^2 + (P(:,2)-1).^2 ));
% x3 = P(index3,1); y3 = P(index3,2);
% P(index3,1) = 0; P(index3,2) = 0;
%
% xmin = min([x1 x2 x3]); xmax = max([x1 x2 x3]);
% dx = 1-xmin;
%
% P(index1,1) = xmin; P(index1,2) = 1+dx;
% P(index2,1) = xmin; P(index2,2) = xmin;
% P(index3,1) = 1+dx; P(index3,2) = xmin;
%
% figure(1);
% plot(x1,y1,'*r'); hold on;
% plot(x2,y2,'*r');
% plot(x3,y3,'*r');
%
% plot(P(index1,1),P(index1,2),'*b'); hold on;
% plot(P(index2,1),P(index2,2),'*b');
% plot(P(index3,1),P(index3,2),'*b');
NElem = help;
xmin = 0.85; dx = 0.15;
P(NElem-2,1) = xmin; P(NElem-2,2) = 1+dx;
P(NElem-1,1) = xmin; P(NElem-1,2) = xmin;
P(NElem,1) = 1+dx; P(NElem,2) = xmin;
figure(1); plot(P(1,1), P(1,2), '*r'); hold on;
for i=1:NElem
    plot(P(i,1), P(i,2), '*r');
end
%

%--------------------------------------------------------- REFLECT POINTSET
function R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha)
eps=1e-8; eta=0.9;
d = Domain('Dist',P);
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+repmat([eps,0],NElem,1))-d)/eps;
n2 = (Domain('Dist',P+repmat([0,eps],NElem,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,2)
R_P(:,1) = P1(I)-2*n1(I).*d(I);
R_P(:,2) = P2(I)-2*n2(I).*d(I);
d_R_P = Domain('Dist',R_P);
J = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
R_P=R_P(J,:); R_P=unique(R_P,'rows');
%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem)
Pc=zeros(NElem,2); A=zeros(NElem,1);
for el = 1:NElem
    vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(Element{el});
    vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
    temp = vx.*vyS - vy.*vxS;
    A(el) = 0.5*sum(temp);
    Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
%------------------------------------------------------- EXTRACT MESH NODES
function [Node,Element] = PolyMshr_ExtrNds(NElem,Node0,Element0)
map = unique([Element0{1:NElem}]);
cNode = 1:size(Node0,1);
cNode(setdiff(cNode,map)) = max(map);
[Node,Element] = PolyMshr_RbldLists(Node0,Element0(1:NElem),cNode);
%----------------------------------------------------- COLLAPSE SMALL EDGES
function [Node0,Element0] = PolyMshr_CllpsEdgs(Node0,Element0,Tol)
while(true)
    cEdge = [];
    for el=1:size(Element0,1)
        if size(Element0{el},2)<4, continue; end;  %Cannot collapse triangles
        vx=Node0(Element0{el},1); vy=Node0(Element0{el},2); nv=length(vx);
        beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
        beta = mod(beta([2:end 1]) -beta,2*pi);
        betaIdeal = 2*pi/size(Element0{el},2);
        Edge = [Element0{el}',Element0{el}([2:end 1])'];
        cEdge = [cEdge; Edge(beta<Tol*betaIdeal,:)];
    end
    if (size(cEdge,1)==0), break; end
    cEdge = unique(sort(cEdge,2),'rows');
    cNode = 1:size(Node0,1);
    for i=1:size(cEdge,1)
        cNode(cEdge(i,2)) = cNode(cEdge(i,1));
    end
    [Node0,Element0] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
%--------------------------------------------------------- RESEQUENSE NODES
function [Node,Element] = PolyMshr_RsqsNds(Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght=cellfun(@length,Element0); nn=sum(ElemLnght.^2);
i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
for el = 1:NElem0
    eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
    i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
    j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
    s(ElemSet) = 1;
    index = index+ElemLnght(el)^2;
end
K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode);
%------------------------------------------------------------ REBUILD LISTS
function [Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode)
Element = cell(size(Element0,1),1);
[foo,ix,jx] = unique(cNode);
if size(Node0,1)>length(ix), ix(end)=max(cNode); end;
Node = Node0(ix,:);
for el=1:size(Element0,1)
    Element{el} = unique(jx(Element0{el}));
    vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(vx);
    [foo,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
    Element{el} = Element{el}(iix);
end
