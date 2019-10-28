function x = projsplx(y)
% Project an n-vector y to the simplex Dn
% Dn = { x : x \in R^n, 0 <= x <= 1, sum(x) = 1}

% Copyright 2019 Daniel Silvestre
% This file is part of OPTool.
%
% OPTool is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% OPTool is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with OPTool.  If not, see <https://www.gnu.org/licenses/>.

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