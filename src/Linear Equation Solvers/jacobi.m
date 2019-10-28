function [ x , parameters] = jacobi( x , ~ , A , b, parameters, ~)
%JACOBI Solves the linear equation A x = b applying the jacobi method 
%   Executes x(k+1) = D^-1 ( b - R x(k) ) 
%   for A = D + R with D a diagonal matrix.

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

% Partition matrix A = D + R
D = diag(diag(A));
R = A - D;

% Residual vector
residual = b - R * x;

% Jacobi iteration
x =  diag(1./diag(D)) * residual;

end

