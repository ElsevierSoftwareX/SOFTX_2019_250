function [ x , parameters] = aor( x , ~ , A , b, parameters, ~)
%AOR Solves the linear equation A x = b applying the Accelerated Over-relaxation method 
%   Executes x(k+1) = (I-rL)^-1 ( (1-w) I + (w-r)L + wU ) x(k) + wc
%   for L = D^-1 A_L, U = D^-1 A_U, c = D^-1 b, and A = D - A_L - A_U 
%   with A_L a strict lower triangular matrix, U a strict upper triangular
%   matrix and D a diagonal matrix.

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

if ~isfield(parameters, 'U') || ~isfield(parameters, 'L')
    % Compute the partition of matrix A
    D = diag(diag(A));
    A_U = -triu(A,1);
    A_L = -tril(A,-1);
    % Compute the U and L matrices
    L = diag(1./diag(D)) * A_L;
    U = diag(1./diag(D)) * A_U;
    % Compute vector c
    c = diag(1./diag(D)) * b;
    % Store U, L and c
    parameters.U = U;
    parameters.L = L;
    parameters.c = c;
else
    % Initialize U, L and c with the values in the parameters
    L = parameters.L;
    U = parameters.U;
    c = parameters.c;
end
% load parameters r and w
r = parameters.r;
w = parameters.w;

% AOR iteration
x = (eye(length(x)) - r * L) \ ( ((1 - w) * eye(length(x)) + (w - r) * L + w * U) * x + w * c );
    
