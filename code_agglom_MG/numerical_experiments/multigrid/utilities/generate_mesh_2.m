function [region]=generate_mesh_2(Dati,step_refinement)
%--------------------------------------------------------------------
% PURPOSE:
%
% This routine generates a structure containing all the information about the
% mesh, i.e.,
%
% region.dim ->  dimension (2)
% region.type_mesh -> type of mesh
% femregion.domain -> (2x2 matrix, real) domain limits
% region.h ->  mesh size
% region.coord -> coordinates of the mesh nodes
% region.connectivity -> connectivity matrix
% region.coords_element -> coordinates of elements (counted with their multiplicity)
% region.boundary_edges -> boundary edges.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

[g, coord, boundary_edge, connectivity]=initialize_mesh(Dati.domain, Dati.type_mesh);

switch Dati.type_mesh
    case{'TRIA_S','TRIA_U'}
        nedge=3;
    case{'QUAD'}
        nedge=4;
end


for i= 1:step_refinement % loop over the number of refinements
    if i==1
        switch Dati.type_mesh
            case{'TRIA_S'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'regular');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'TRIA_U'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'longest');
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'longest');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'QUAD'}
                [coord, connectivity]= refine_quad(coord, connectivity);
        end
        
    else
        switch Dati.type_mesh
            case{'TRIA_S','TRIA_U'}
                [coord, boundary_edge, connectivity]=refine_mesh(g,coord, boundary_edge, connectivity,'regular');
                coord=jiggle_mesh(coord, boundary_edge, connectivity);
            case{'QUAD'}
                [coord, connectivity]= refine_quad(coord, connectivity);
        end
    end
end
ne=size(connectivity,2);
v = [Dati.domain(1,1) Dati.domain(2,1);Dati.domain(1,2) Dati.domain(2,1);Dati.domain(1,2) Dati.domain(2,2);Dati.domain(1,1) Dati.domain(2,2)];
h= sqrt(polyarea(v(:,1),v(:,2))/ne);

connectivity = connectivity(1:nedge,:);

coords_element = []; % coordinates of the elements (counted with their multiplicity)
%Build BBox (add by GP)
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
    'type_mesh',Dati.type_mesh,...
    'domain',Dati.domain,...
    'nedges', nedges,...
    'h',h,...
    'nvert',size(coord,2),...
    'ne',ne,...
    'coord',coord',...
    'boundary_edges',boundary_edge,...
    'connectivity',connectivity,...
    'coords_element',coords_element);

%Added by GP
region.BBox = BBox;
region.area = elem_area;
region.max_kb = max_kb;
region.connectivity = Element;
region.coords_element = coords_element_cell;

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
%region.h = max(H);
region.h = H;

end


function [g, coord, boundary_edge, connectivity]=initialize_mesh(domain, type_mesh)
%--------------------------------------------------------------------
% PURPOSE:
%
% This routine initializes a tringular mesh on a given domain
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

% domain
x0=domain(1,1);
x1=domain(1,2);
y0=domain(2,1);
y1=domain(2,2);

% geometry
g=[
    2     2     2     2
    x0    x1    x1    x0
    x1    x1    x0    x0
    y1    y1    y0    y0
    y1    y0    y0    y1
    0     0     0     0
    1     1     1     1
    ];

% points
coord=[
    x0     x1     x1    x0
    y1     y1     y0    y0
    ];

% boundary edges
boundary_edge=[
    1     2     3     4
    2     3     4     1
    0     0     0     0
    1     1     1     1
    1     2     3     4
    0     0     0     0
    1     1     1     1
    ];

% connectivity
switch type_mesh
    case{'TRIA_S','TRIA_U'}
        connectivity =[
            1     2
            4     1
            3     3
            1     1
            ];
    case{'QUAD'}
        connectivity =[
            1
            4
            3
            2
            1
            ];
end

end




function [p1,e1,t1,u1] =  refine_mesh(g,p,e,t,u,it,mode)
%--------------------------------------------------------------------
% PURPOSE:
%
% This routine refines a triangular mesh.
% [p1,e1,t1]=refine_mesh(g,p,e,t) returns a refined version
% of the triangular mesh specified by the geometry g, point matrix p,
% edge matrix e, and triangle matrix t.
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

np=size(p,2);
ne=size(e,2);
nt=size(t,2);

if nargout==4
    intp=1;
else
    intp=0;
end

if nargin-intp==4
    it=(1:nt)';                           % All triangles
    mode='regular';
end

if (~intp) && nargin==5
    it=u;
    if ischar(it)
        mode=it;
        it=(1:nt)';                         % All triangles
    else
        mode='regular';
    end
end

if (~intp) && nargin==6
    mode=it;
    it=u;
end

if intp && nargin==6
    if ischar(it)
        mode=it;
        it=(1:nt)';                         % All triangles
    else
        mode='regular';
    end
end

if strcmp(mode,'regular')==0 && strcmp(mode,'longest')==0
    error('PDE:refinemesh:InvalidRefineMode', 'Unknown refinement mode.');
end

if size(it,1)>1                        % Triangles
    it=it';
else                                    % Subdomains
    it=pdesdt(t,it);
end

% Cannot use matrix indices that exceeds the size of a signed int
[comp,maxsize]=computer;
indexproblem=np^2>maxsize;

% Find longest side of each triangle
ls=3*ones(1,nt);
d1=(p(1,t(1,:))-p(1,t(2,:))).^2+(p(2,t(1,:))-p(2,t(2,:))).^2;
d=(p(1,t(2,:))-p(1,t(3,:))).^2+(p(2,t(2,:))-p(2,t(3,:))).^2;
ii=find(d>d1);
ls(ii)=1*ones(size(ii));
d1=max(d,d1);
d=(p(1,t(3,:))-p(1,t(1,:))).^2+(p(2,t(3,:))-p(2,t(1,:))).^2;
ii=find(d>d1);
ls(ii)=2*ones(size(ii));
% Permute so longest side is 3
ii=find(ls==1);
d=t(1,ii);
t(1,ii)=t(2,ii);
t(2,ii)=t(3,ii);
t(3,ii)=d;
ii=find(ls==2);
d=t(1,ii);
t(1,ii)=t(3,ii);
t(3,ii)=t(2,ii);
t(2,ii)=d;

itt1=ones(1,nt);
itt1(it)=zeros(size(it));
it1=find(itt1);                         % Triangles not yet to be refined
it=find(itt1==0);                       % Triangles whos longest side is to be bisected

% Make a connectivity matrix, with edges to be refined.
% -1 means no point is yet allocated
ip1=t(1,it);
ip2=t(2,it);
if strcmp(mode,'regular')
    ip3=t(3,it);
end
A=sparse(ip1,ip2,-1,np,np);
if strcmp(mode,'regular')
    A=A+sparse(ip2,ip3,-1,np,np);
    A=A+sparse(ip3,ip1,-1,np,np);
end
A=-((A+A.')<0);
newpoints=1;

% loop until no additional hanging nodes are introduced
while newpoints
    newpoints=0;
    n=length(it1);
    ip1=t(1,it1);
    ip2=t(2,it1);
    ip3=t(3,it1);
    m1=zeros(1,n);
    m2=m1;
    m3=m1;
    for i=1:n
        m3(i)=A(ip1(i),ip2(i));
        m1(i)=A(ip2(i),ip3(i));
        m2(i)=A(ip3(i),ip1(i));
    end
    ii=find(m3);
    if length(ii)>0
        itt1(it1(ii))=zeros(size(ii));
    end
    ii=find((m1 | m2) & (~m3));
    if length(ii)>0
        A=A+sparse(ip1(ii),ip2(ii),-1,np,np);
        A=-((A+A.')<0);
        newpoints=1;
        itt1(it1(ii))=zeros(size(ii));
    end
    it1=find(itt1);                       % Triangles not yet fully refined
    it=find(itt1==0);                     % Triangles fully refined
end

% Find edges to be refined
if ~indexproblem
    ie=full(A(e(1,:)+(e(2,:)-1)*np))==-1;
else
    ie=l_extract(A,e(1,:),e(2,:))==-1;
end

ie1=find(ie==0);                        % Edges not to be refined
ie=find(ie);                            % Edges to be refined

% Get the edge "midpoint" coordinates
[x,y]=pdeigeom(g,e(5,ie),(e(3,ie)+e(4,ie))/2);
% Create new points
p1=[p [x;y]];
if intp
    u1=[u;(u(e(1,ie),:)+u(e(2,ie),:))/2];
end
ip=(np+1):(np+length(ie));
np1=np+length(ie);
% Create new edges
e1=[e(:,ie1) ...
    [e(1,ie);ip;e(3,ie);(e(3,ie)+e(4,ie))/2;e(5:7,ie)] ...
    [ip;e(2,ie);(e(3,ie)+e(4,ie))/2;e(4,ie);e(5:7,ie)]];
% Fill in the new points
if ~indexproblem
    A(e(1,ie)+np*(e(2,ie)-1))=ip;
    A(e(2,ie)+np*(e(1,ie)-1))=ip;
else
    A=l_assign(A,[e(1,ie) e(2,ie)],[e(2,ie) e(1,ie)],[ip ip]);
end

% Generate points on interior edges
[i1,i2]=find(A==-1 & A.'==-1);
i=find(i2>i1);
i1=i1(i);
i2=i2(i);
p1=[p1 ((p(1:2,i1)+p(1:2,i2))/2)];
if intp
    u1=[u1;(u(i1,:)+u(i2,:))/2];
end
ip=(np1+1):(np1+length(i));
np1=np1+length(i);
% Fill in the new points
if ~indexproblem
    A(i1+np*(i2-1))=ip;
    A(i2+np*(i1-1))=ip;
else
    A=l_assign(A,[i1 i2],[i2 i1],[ip ip]);
end

% Lastly form the triangles
ip1=t(1,it);
ip2=t(2,it);
ip3=t(3,it);
if ~indexproblem
    mp1=full(A(ip2+np*(ip3-1)));
    mp2=full(A(ip3+np*(ip1-1)));
    mp3=full(A(ip1+np*(ip2-1)));
else
    mp1=l_extract(A,ip2,ip3);
    mp2=l_extract(A,ip3,ip1);
    mp3=l_extract(A,ip1,ip2);
end

% Find out which sides are refined
bm=1*(mp1>0)+2*(mp2>0);
% The number of new triangles
nt1=length(it1)+length(it)+sum(mp1>0)+sum(mp2>0)+sum(mp3>0);
t1=zeros(4,nt1);
t1(:,1:length(it1))=t(:,it1);           % The unrefined triangles
nnt1=length(it1);
if isempty(bm)
    i = bm;
else
    i=find(bm==3);                          % All sides are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[mp1(i);mp2(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
    i = bm;
else
    i=find(bm==2);                          % Sides 2 and 3 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
    i = bm;
else
    i=find(bm==1);                          % Sides 3 and 1 are refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp3(i);mp1(i);t(4,it(i))];
nnt1=nnt1+length(i);
if isempty(bm)
    i = bm;
else
    i=find(bm==0);                          % Side 3 is refined
end
t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
nnt1=nnt1+length(i);
end % end of function refine_mesh
%--------------------------------------------------------------------
%--------------------------------------------------------------------


function k=l_extract(A,i,j)

if numel(i)~=numel(j)
    error('PDE:refinemesh:ijNumel', 'i and j must have the same number of elements.')
end

k=zeros(size(i));

for l=1:numel(i)
    k(l)=A(i(l),j(l));
end
end
%--------------------------------------------------------------------
%--------------------------------------------------------------------


function A=l_assign(A,i,j,k)

if numel(i)~=numel(j) || numel(i)~=numel(k)
    error('PDE:refinemesh:ijkNumel', 'i, j, and k must have the same number of elements.')
end

for l=1:numel(i)
    A(i(l),j(l))=k(l);
end

end

%--------------------------------------------------------------------
%--------------------------------------------------------------------
function [x,y]=pdeigeom(dl,bs,s)

if nargin<1 || nargin>3
    error('PDE:pdeigeom:nargin', 'pdeigeom should have 1-3 input arguments.');
end

if ischar(dl)
    if nargin==1
        x=feval(dl);
    elseif nargin==2
        x=feval(dl,bs);
    else                                  % nargin==3
        [x,y]=feval(dl,bs,s);
    end
    return
end

nbs=size(dl,2);

if nargin==1
    x=nbs;                                % number of boundary segments
    return
end

d=[zeros(1,nbs);
    ones(1,nbs);
    dl(6:7,:)];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs)
    error('PDE:pdeigeom:InvalidBs', 'Non-existent boundary segment number.')
end

if nargin==2
    x=d(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1
    bs=bs*ones(size(s));                  % expand bs
elseif m~=size(s,1) || n~=size(s,2)
    error('PDE:pdeigeom:SizeBs', 'bs must be scalar or of same size as s.');
end

if ~isempty(s)
    for k=1:nbs
        ii=find(bs==k);
        if ~isempty(ii)
            x0=dl(2,k);
            x1=dl(3,k);
            y0=dl(4,k);
            y1=dl(5,k);
            if dl(1,k)==1                     % Circle fragment
                xc=dl(8,k);
                yc=dl(9,k);
                r=dl(10,k);
                a0=atan2(y0-yc,x0-xc);
                a1=atan2(y1-yc,x1-xc);
                if a0>a1
                    a0=a0-2*pi;
                end
                theta=(a1-a0)*(s(ii)-d(1,k))/(d(2,k)-d(1,k))+a0;
                x(ii)=r*cos(theta)+xc;
                y(ii)=r*sin(theta)+yc;
            elseif dl(1,k)==2                 % Line fragment
                x(ii)=(x1-x0)*(s(ii)-d(1,k))/(d(2,k)-d(1,k))+x0;
                y(ii)=(y1-y0)*(s(ii)-d(1,k))/(d(2,k)-d(1,k))+y0;
            elseif dl(1,k)==4                 % Ellipse fragment
                xc=dl(8,k);
                yc=dl(9,k);
                r1=dl(10,k);
                r2=dl(11,k);
                phi=dl(12,k);
                t=[r1*cos(phi) -r2*sin(phi); r1*sin(phi) r2*cos(phi)];
                rr0=t\[x0-xc;y0-yc];
                a0=atan2(rr0(2),rr0(1));
                rr1=t\[x1-xc;y1-yc];
                a1=atan2(rr1(2),rr1(1));
                if a0>a1
                    a0=a0-2*pi;
                end
                % s should be proportional to arc length
                % Numerical integration and linear interpolation is used
                nth=100;                % The number of points in the interpolation
                th=linspace(a0,a1,nth);
                rr=t*[cos(th);sin(th)];
                theta=pdearcl(th,rr,s(ii),d(1,k),d(2,k));
                rr=t*[cos(theta);sin(theta)];
                x(ii)=rr(1,:)+xc;
                y(ii)=rr(2,:)+yc;
            else
                error('PDE:pdeigeom:InvalidSegType', 'Unknown segment type.');
            end
        end
    end
end
end



function p=jiggle_mesh(p,e,t,p1,v1,p2,v2)
%--------------------------------------------------------------------
% PURPOSE:
%
% This routines jiggle internal points of a triangular mesh by adjusting
% the node point positions. The quality of the mesh normally increases.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

% Error checks
nargs = nargin;
if rem(nargs+3,2)
    error('PDE:jigglemesh:NoParamPairs','Param value pairs expected.')
end

% Default values
Opt='off';
Iter=-1;

for i=4:2:nargs
    Param = eval(['p' int2str((i-4)/2 +1)]);
    Value = eval(['v' int2str((i-4)/2 +1)]);
    if ~ischar(Param)
        error('PDE:jigglemesh:ParamNotString', 'Parameter must be a string.')
    elseif size(Param,1)~=1
        error('PDE:jigglemesh:ParamEmptyOrNot1row', 'Parameter must be a non-empty single row string.')
    end
    Param = lower(Param);
    if strcmp(Param,'opt')
        Opt=lower(Value);
        if ~ischar(Opt)
            error('PDE:jigglemesh:OptNotString', 'Opt must be a string.')
        elseif ~strcmp(Opt,'off') && ~strcmp(Opt,'minimum') && ~strcmp(Opt,'mean')
            error('PDE:jigglemesh:OptInvalidString', 'Opt must be off | minimum | {mean}.')
        end
    elseif strcmp(Param,'iter')
        Iter=Value;
        if ischar(Iter)
            error('PDE:jigglemesh:IterString', 'Iter must not be a string.')
        elseif ~all(size(Iter)==[1 1])
            error('PDE:jigglemesh:IterNotScalar', 'Iter must be a scalar.')
        elseif imag(Iter)
            error('PDE:jigglemesh:IterComplex', 'Iter must not be complex.')
        elseif Iter<-1
            error('PDE:jigglemesh:IterNeg', 'Iter must be non negative.')
        end
    else
        error('PDE:jigglemesh:InvalidParam', ['Unknown parameter: ' Param])
    end
end

if Iter==-1 && strcmp(Opt,'off')
    Iter=1;
elseif Iter==-1
    Iter=20;
end

% Determine interior non boundary points
ep=sort([e(1,:) e(2,:)]);
i=ones(1,size(p,2));
j=ep(find([1 sign(diff(ep))]));
i(j)=zeros(size(j));
i=find(i);

np=size(p,2);
nt=size(t,2);

if ~strcmp(Opt,'off')
    q=pdetriq(p,t);
    if strcmp(Opt,'minimum')
        q=min(q);
    else
        q=mean(q);
    end
end

j=1;
while j<=Iter
    X=sparse(t([1 2 3],:),t([2 3 1],:),p(1,t(1:3,:)),np,np);
    Y=sparse(t([1 2 3],:),t([2 3 1],:),p(2,t(1:3,:)),np,np);
    N=sparse(t([1 2 3],:),t([2 3 1],:),1,np,np);
    m=sum(N);
    X=sum(X)./m;
    Y=sum(Y)./m;
    p1=p;
    p(1,i)=X(i);
    p(2,i)=Y(i);
    if ~strcmp(Opt,'off')
        q1=q;
        q=pdetriq(p,t);
        if strcmp(Opt,'minimum')
            q=min(q);
        elseif strcmp(Opt,'mean')
            q=mean(q);
        end
        if q<q1
            p=p1;
            break,
        elseif q1+1e-4>q
            break
        end
    end
    j=j+1;
end
end


function [pn, tn]= refine_quad(p, t)

xnod = p(1,:);
ynod = p(2,:);

nodes = t;

nele=length(nodes(1,:));
nno=length(xnod);
nodesn=[];
bnodn=[];
sum=ones(1,4)/4.;
check=sparse(1,1,1,nno,nno,10*nno);
check(1,1)=0;

for iel=1:nele
    
    %     new nodes
    iv=nodes(1:end-1,iel);
    xnod=[xnod, sum*xnod(iv)'];
    ynod=[ynod, sum*ynod(iv)'];
    nno=nno+1; nodmid=nno;
    if check(iv(1),iv(2))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(1))+xnod(iv(2)))/2.];
        ynod=[ynod, (ynod(iv(1))+ynod(iv(2)))/2.];
        check(iv(1),iv(2))=nno;check(iv(2),iv(1))=nno;
    end
    
    if check(iv(3),iv(2))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(3))+xnod(iv(2)))/2.];
        ynod=[ynod, (ynod(iv(3))+ynod(iv(2)))/2.];
        check(iv(3),iv(2))=nno;check(iv(2),iv(3))=nno;
    end
    
    if check(iv(3),iv(4))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(3))+xnod(iv(4)))/2.];
        ynod=[ynod, (ynod(iv(3))+ynod(iv(4)))/2.];
        check(iv(3),iv(4))=nno;check(iv(4),iv(3))=nno;
    end
    if check(iv(1),iv(4))==0
        nno=nno+1;
        xnod=[xnod, (xnod(iv(1))+xnod(iv(4)))/2.];
        ynod=[ynod, (ynod(iv(1))+ynod(iv(4)))/2.];
        check(iv(1),iv(4))=nno;check(iv(4),iv(1))=nno;
    end
    ivn=[iv(1) check(iv(1),iv(2)) nodmid check(iv(1),iv(4))];
    nodesn=[nodesn, ivn'];
    ivn=[check(iv(1),iv(2)) iv(2) check(iv(2),iv(3)) nodmid];
    nodesn=[nodesn, ivn'];
    ivn=[nodmid check(iv(3),iv(2)) iv(3) check(iv(3),iv(4))];
    nodesn=[nodesn, ivn'];
    ivn=[check(iv(1),iv(4)) nodmid check(iv(3),iv(4)) iv(4)];
    nodesn=[nodesn, ivn'];
end

pn(1,:)=xnod; pn(2,:)=ynod;

tn = [nodesn; ones(1,length(nodesn))];

end
