function x = projsplx(y)
% Project an n-vector y to the simplex Dn
% Dn = { x : x \in R^n, 0 <= x <= 1, sum(x) = 1}

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for i = 1:m-1
    tmpsum = tmpsum + s(i);
    tmax = (tmpsum - 1)/i;
    if tmax >= s(i+1)
        bget = true;
        break;
    end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end

x = max(y-tmax,0);

return;