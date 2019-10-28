function T = Taor(A,~,w,r)
%TAOR Computes the transition matrix Taor for the AOR algorithm.
%   Computes :
%   T = (I - r * L)^-1 ((1 - w) * I + (w - r) * L + w * U);

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

% Compute the partition of matrix A
D = diag(diag(A));
A_U = -triu(A,1);
A_L = -tril(A,-1);
% Compute the U and L matrices
L = diag(1./diag(D)) * A_L;
U = diag(1./diag(D)) * A_U;

T = (eye(length(A)) - r * L) \ ((1 - w) * eye(length(A)) + (w - r) * L + w * U);

end

