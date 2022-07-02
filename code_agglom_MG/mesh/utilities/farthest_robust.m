function [P,scores] = farthest_robust(vert,k)
if k == 2
    D = pdist(vert);
    M = max(D);
    Z = squareform(D);
    [i1,i2] = find(Z == M,1,'first');
    P = vert([i1,i2],:);
    return
end

% I = randperm(size(vert,1));
% P = vert(I(1:k),:);
% left = vert(I(k+1:end),:);
[P,left] = farthest_first(vert,k);

scores = zeros(1,k);

go_on = true;
while go_on
    go_on = false;
    for i = 1:size(P,1)
        scores(i) = score_fun(P,i);
        jmax = 0;
        for j = 1:size(left,1)
             P_tmp = P;
             P_tmp(i,:) = left(j,:);
             sj = score_fun(P_tmp,i);
             if sj > scores(i)
                 scores(i) = sj;
                 jmax = j;
                 break
             end
        end
        if jmax > 0
            tmp = left(jmax,:);
            left(jmax,:) = P(i,:);
            P(i,:) = tmp;
            go_on = true;
        end
    end
end


end

function s = score_fun(P,i)
try
    [~,s_in] = convhull(P);
catch
    s = 0;
    return
end
P(i,:) = [];
try
    [~,s_out] = convhull(P);
catch
    s = s_in;
    return
end
s = s_in - s_out;
end
