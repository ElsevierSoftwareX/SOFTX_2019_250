function [w_max, r_max] = optAOR(A)
%OPTSOR Returns the optimal w that guarantees the fastest convergence for
%the Successive Over-relaxation method.
%   Solves the optimization problem:
%       min    max{ lambda( T_SOR ) }
%   Given the format of the problem, the solution is found through brute
%   force.

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

% ==== Could NOT solve using optimization (to research) ====
% ==== Brute Force Option ====
% Compute the partition of matrix A
D = diag(diag(A));
A_U = -triu(A,1);
A_L = -tril(A,-1);
% Compute the U and L matrices
L = diag(1./diag(D)) * A_L;
U = diag(1./diag(D)) * A_U;

N = 100;
w_values = linspace(0,2,N);
r_values = linspace(0,2,N);
spectral_radius = zeros(N,N);

% In case it is not possible select the Jacobi parameters
w_max = 1;
r_max = 0;

rho_min = 1;

for i = 1:N
    for j = 1:N
        w = w_values(i);
        r = r_values(j);
        T = (eye(length(A)) - r * L) \ ((1 - w) * eye(length(A)) + (w - r) * L + w * U);
        
        lambdaT = eig(T);
        
        spectral_radius(i,j) = max(abs(lambdaT));
        
        if spectral_radius(i,j) < rho_min
            rho_min = spectral_radius(i,j);
            w_max = w;
            r_max = r;
        end
    end
end

% surf(w_values,r_values,spectral_radius);

% ==== Closed-form solution ====
%


