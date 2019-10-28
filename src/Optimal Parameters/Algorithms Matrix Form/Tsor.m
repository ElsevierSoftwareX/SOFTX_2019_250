function T = Tsor(A,~,w)
%TSOR Computes the transition matrix Tsor for the SOR algorithm.
%   Computes :
%   T = -(D+(w*L))^-1*(w*U + (w-1)*D) 
% where we explore the fact that:
%   iLw = (D+(w*L))^-1 = \sum_{j=0}^{n-1} (-1)^j D*(w*L)^j

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

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

% iLw = D;
% for j = 1:length(A)-1
%     iLw = iLw + (-1)^j* D*(w*L)^j;
% end
% T = -iLw*(w*U + (w-1)*D);

T = -inv(D+(w*L))*(w*U + (w-1)*D);

end

