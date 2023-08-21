%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Modification of PolyMesher      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [region] = generate_mesh_L_shaped(Domain,NElem,MaxIter,P)
flagFixedSeeds = 0;
if ~exist('P','var')
    %P=PolyMshr_RndPtSet(NElem,Domain);
    %P=PolyMshr_NormalRndPtSet(NElem,Domain);
    P=PolyMshr_FixedPtSet(NElem,Domain); flagFixedSeeds = 1;
end
NElem = size(P,1);
Tol=5e-3; It=0; Err=1; c=1.5;%c=1.5;
out = Domain('BdBoxes'); BdBoxA = out{1}; BdBoxB = out{2};
%Area = (BdBoxA(2)-BdBoxA(1))*(BdBoxA(4)-BdBoxA(3)) - (BdBoxB(2)-BdBoxB(1))*(BdBoxB(4)-BdBoxB(3));
Area = (BdBoxA(2)-BdBoxA(1))*(BdBoxA(4)-BdBoxA(3));
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

if(flagFixedSeeds)
    
    % delete node 37, and modify node 41 as (1.0, 1.0)
    Element{6}(7) = [];
    Element{10}(4) = 41;
    Node(41,1) = 1; Node(41,2) = 1;
    
    %delete node 37
    Node(37,:) = [];
    ne = length(Element);
    for i=1:ne
        for j=1:length(Element{i})
            if(Element{i}(j)>37)
                Element{i}(j) = Element{i}(j) - 1;
            end
        end
    end
    
    
    % there are some nodes very close each other, here we delete them from
    % connectivity
    toll = 1e-6;
    modify = []; count = 0;
    Nnodes = size(Node,1);
    for i = 2:size(Node,1)
        if( sqrt( (Node(i,1)-Node(i-1,1))^2 + (Node(i,2)-Node(i-1,2))^2) < toll )
            count = count + 1;
            modify = [modify, i];
        end
    end
    
    %Node(modify,:) = [];
    for k=modify
        for ie=1:ne
            I = find(Element{ie}==k);
            Element{ie}(I) = k-1;
        end
    end
    
    % queste coppie di nodi ( (12,38) e (2,9) ) sono uguali ma non
    % consecutivi: li correggo a mano perché non c'ho più sbatti!
    for ie = 1:ne
        for j=1:length(Element{ie})
            if(Element{ie}(j)==38)
                Element{ie}(j)=12;
            end
        end
    end
    
    for ie = 1:ne
        for j=1:length(Element{ie})
            if(Element{ie}(j)==9)
                Element{ie}(j)=2;
            end
        end
    end
    
    for ie=1:ne
        delete = [];
        for j=2:length(Element{ie})
            if(Element{ie}(j)==Element{ie}(j-1))
                delete = [delete, j];
            end
        end
        Element{ie}(delete) = [];
    end
    
    Node(21,1) = Node(21,1) - 0.15;
    Node(37,1) = Node(37,1) - 0.15;
    Node(37,2) = Node(37,2) + 0.01;
    Node(35,2) = Node(35,2) + 0.1;
    
end


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

region.tag = [1:region.ne];

%per capire quali sono le coordinate degli elementi dell'elemento
%coords_element{length(Element)}

%%% POLYMESHER FUNCTIONS %%%%

%------------------------------------------------- GENERATE RANDOM POINTSET
function P = PolyMshr_RndPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem
    Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
    Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
    d = Domain('Dist',Y);
    I = find(d(:,end)<0);                 %Index of seeds inside the domain
    NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
    P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
    Ctr = Ctr+NumAdded;
end
%figure(2);
%plot(P(:,1),P(:,2),'*'); axis equal;

% this generate point with non uniform distributin, but it does not work
% because you should define also a size function h(x) somewhere
function P = PolyMshr_NormalRndPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem
    
    pd = makedist('Triangular','a',0,'b',0.5,'c',1);
    for k=1:NElem
        x(k) = random(pd);
        y(k) = random(pd);
    end
    
    %x = normrnd(0,1,[NElem,1]);
    %xmin = min(x); xmax = max(x);
    %x = (1/(xmax-xmin))*x - xmin/(xmax-xmin);
    Y(:,1) = (BdBox(2)-BdBox(1))*x+BdBox(1);
    
    %y = normrnd(0,1,[NElem,1]);
    %ymin = min(y); ymax = max(y);
    %y = (1/(ymax-ymin))*y - ymin/(ymax-ymin);
    Y(:,2) = (BdBox(4)-BdBox(3))*y+BdBox(3);
    d = Domain('Dist',Y);
    I = find(d(:,end)<0);                 %Index of seeds inside the domain
    NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
    P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
    Ctr = Ctr+NumAdded;
end
figure(2);
plot(P(:,1),P(:,2),'*'); axis equal;

function P = PolyMshr_FixedPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
dx = (BdBox(2)-BdBox(1))/8; eps = dx/3;
% P = [dx dx
%     3*dx+eps   dx
%     5*dx       dx+2*eps
%     7*dx-eps   dx-eps
%     dx+eps     3*dx
%     3*dx       3*dx % fixed
%     5*dx       3*dx % fixed
%     7*dx       3*dx
%     dx+2*eps   5*dx
%     3*dx       5*dx % fixed
%     dx         7*dx+eps
%     3*dx       7*dx
%      2*dx       4*dx-eps
%      dx         6*dx
%      5*dx       dx-eps
%      2*dx-eps   dx-eps
%     ];


P = [dx dx
    3*dx+eps   dx
    5*dx       dx+2*eps
    7*dx-eps   dx-eps
    dx+eps     3*dx
    3*dx       3*dx % fixed
    5*dx       3*dx % fixed
    7*dx       3*dx
    dx+2*eps   5*dx
    3*dx       5*dx % fixed
    dx         7*dx+eps
    3*dx       7*dx
    2*dx       4*dx-eps
    dx         6*dx
    5*dx       dx-eps
    2*dx-eps   dx-eps
    ];
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
