function [P,left] = farthest_first(vert,k)
if k == 2
    D = pdist(vert);
    M = max(D);
    Z = squareform(D);
    [i1,i2] = find(Z == M,1,'first');
    P = vert([i1,i2],:);
    return
end


% check wikipedia farthest-first traversal

P = zeros(k+1,size(vert,2));
try
    P(1,:) = conv_centr(vert);
catch
    P(1,:) = mean(vert,1);
end
left = vert;
N = size(left,1);
% k = min(k,N);
for i = 1:k
    d = 0;
    for j = 1:(N-i+1)
        dnew = prod(vecnorm(P(1:i,:)-left(j,:),2,2));
        if dnew > d
            d = dnew;
            imax = j;
        end
    end
    P(i+1,:) = left(imax,:);
    left(imax,:) = [];
end
P(1,:) = [];
end