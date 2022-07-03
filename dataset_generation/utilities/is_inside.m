function [TFin,TFon] = is_inside(view1,points,convex)
if nargin < 3
    convex = view1.convex;
end

if manifold(view1) == 3
    vert = vertices(view1);
    shp = alphaShape(vert,inf);
    TFconv = inShape(shp,points(:,1),points(:,2),points(:,3));
    if convex
        TFin = TFconv;
        TFon = [];
        return
    else
        N = size(points,1);
        points = points(TFconv,:);
        [TFinconv,TFonconv] = inside_fun(view1,points);
        TFin = false(1,N);
        TFon = false(1,N);
        TFin(TFconv)=TFinconv;
        TFon(TFconv)=TFonconv;
        return
    end
end
[TFin,TFon] = inside_fun(view1,points);
end

function [TFin,TFon] = inside_fun(view1,points)
% view1 must contain only ONE element

dim = manifold(view1);
if dim == 0
    TFin = vecnorm(points - get0(view1,1),inf,2) < view1.tol;   
    TFon = TFin;
    return
end

N = size(points,1);
done = false(1,N);
TFin = false(1,N);
TFon = false(1,N);
toldim = view1.tol*view1.dim;
vert =  vertices(view1);
v0 = vert -  vert(1,:);
p0 = points - vert(1,:);
rv0 = rank(v0,toldim);
for j = 1:N
    % if point is not on element plane
    done(j) =  rv0 < rank([v0;p0(j,:)],toldim);
end

if all(done)
    return
end

I = find(not(done));
N = length(I);
% find random direction from point inside the element plane
w = rand(1,view1.elem_num(1));
v = w/sum(w)*v0;
% should check that v is not zero and not parallel to any element face
% but these events have 0 probability

t_list = zeros(N,view1.elem_num(dim)); % line = t*v+point

for i = 1:view1.elem_num(dim) % for each element face
    face = RView(view1,i,dim-1);
    vert = vertices(face);

    % face on the new origin --> face plane = vertices combination
    p0 = points(I,:) - vert(1,:);
    v0 = vert(2:end,:) - vert(1,:);

    % v-line from p0 : t*v + p0 = v0'*x : face plane
    % p0 = [vert;v]'*[x;-t] = A*y

    q = zeros(N,view1.dim);
    for j = 1:N
        y = [v0;v]'\p0(j,:)';
        t_list(j,i) = -y(end);
        q(j,:) = t_list(j,i)*v + points(I(j),:); % intersection line-plane
    end

    t_list(not(inside_fun(face,q)),i) = nan;
    uon = abs(t_list(:,i)) < view1.tol;
    TFon(I(uon)) = true;
    TFin(I(uon)) = true;
    done(I(uon)) = true;
end

if all(done)
    return
end

for j = find(not(done(I)))
    % using the list of intersections along the line
    % find if the point is inside the element
    t = t_list(j,:);
    t(isnan(t)) = [];
    if not(isempty(t_list(j,:)))
        t = uniquetol(t,view1.tol);
        ID = find(t > 0, 1, 'first');        
        TFin(I(j)) = not(isempty(ID)) && (mod(ID,2)==0);
    end
end

end