function [ x , parameters] = quasiChebyshev( x , previous_x , A , b, parameters, ~)
%QUASICHEBYSHEV Solves the linear equation A x = b applying the Quasi-Chebyshev accelaration
%for the Jacobi method as given in [4]
%   Executes x(k+1) = w_(k+1) *(x_jacobi(k+1) - x(k-1) ) + x(k-1)
%   For details on the computation of w_(k+1) see [4] [Wen] Quasi-Chebyshev
%   accelerated iteration methods based on optimization for linear systems.

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

if ~isfield(parameters, 'M') || ~isfield(parameters, 'N')
    if parameters.jacobi == 1
        M = diag(diag(A));
        N = M - A;
    else
        % Parameter alpha > 0
        alpha = parameters.alpha;
        
        % Partition matrix A = M - N using H and S
        H = 0.5 * (A + A');
        S = 0.5 * (A - A');
        M = 1/(2*alpha) * (alpha * eye(length(x)) + H) * (alpha * eye(length(x)) + S);
        N = 1/(2*alpha) * (alpha * eye(length(x)) - H) * (alpha * eye(length(x)) - S);
    end
    % Store all variables for the next iteration
    parameters.M = M;
    parameters.N = N;
end

% If the residual is too small causes problems with the definition of w
if norm(b-A*x) < eps
    return
end

% Load variables
M = parameters.M;
N = parameters.N;

% Compute x_jacobi(k+1)
x_jacobi = diag(1./diag(M)) * (N * x + b);

[~, p] = chol(A);

if p == 0 && norm(A-A') < 1e-5
    % A is symmetric positive definite
    w = ( (x_jacobi - previous_x)' * A * (x_jacobi - previous_x) )^-1 * (b - A * previous_x)' * (x_jacobi - previous_x);
else
    % A is not symmetric positive definite
    w = ( (A * x_jacobi - b -(A * previous_x -b))' * (eye(length(x)) + 0.5 * (A + A'))^-2 * (A * x_jacobi - b -(A * previous_x -b)))^-1 * (b - A * previous_x)' * (eye(length(x)) + 0.5 * (A + A'))^-2 * (A * x_jacobi - b -(A * previous_x -b));
end

% Quasi-Chebyshev iteration
x =  w * (x_jacobi  - previous_x) + previous_x;

end

