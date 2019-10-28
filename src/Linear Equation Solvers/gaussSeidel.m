function [ x , parameters] = gaussSeidel( x , ~ , A , b, parameters, ~)
%gaussSeidel Solves the linear equation A x = b applying the Gauss-Seidel method 
%   Executes x(k+1) = L^-1 ( b - U x(k) ) 
%   for A = L + U with L a lower triangular matrix and U a strict upper
%   triangular matrix. To avoid the inverse, the element-wise definition is
%   used: 
%   x_i(k+1) = 1/a_ii (b_i - \sum_{j=1}^{i-1} a_ij x_j(k+1) - \sum_{j=i+1}^{n} a_ij x_j(k))

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

n = length(x);

for i = 1:n
    % Gauss-Seidel element-wise iteration
    x(i) = (1/A(i, i))*( b(i) - A(i, [1:i-1 i+1:end])*x([1:i-1 i+1:end],1) );
end

