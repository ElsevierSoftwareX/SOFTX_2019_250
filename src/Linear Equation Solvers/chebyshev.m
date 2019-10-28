function [ x , parameters] = chebyshev( x , previous_x , A , b, parameters, ~)
%Chebyshev Solves the linear equation A x = b applying the Chebyshev accelaration
%for the Jacobi method as given in [4]
%   Executes x(k+1) = w_(k+1) *(T x(k) + M^-1 b - x(k) ) + x(k) 
%   for T = M^-1 N and A = M - N. 
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

if ~isfield(parameters, 'rho')
    % Partition matrix A = M - N
    M = diag(diag(A));
    N = M - A;
    
    % Compute the transition matrix
    T = M \ N;
    % Store the spectral radius of T
    rho = max(abs(eig(T)));
    
    % Compute the first w
    w = 2/(2 - rho^2);
    
    % Compute the first auxiliary variables C1 and C2
    Cm = [1/rho 2/rho^2-1];
    
    % Store all variables for the next iteration
    parameters.rho = rho;
    parameters.Cm = Cm;
    parameters.T = T;
    parameters.M = M;
    
else
    % Load variables
    rho = parameters.rho;
    Cm = parameters.Cm;
    T = parameters.T;
    M = parameters.M;
    
    % Compute the next Cm using the recurrence formula
    C_next_m = 2/rho * Cm(2) - Cm(1);
    
    % Next parameter w
    w = (2 * Cm(2) ) / (rho * C_next_m);
    
    % Save the auxiliary variable Cm
    Cm = [Cm(2) C_next_m];
    parameters.Cm = Cm;
end
    
% Chebyshev iteration
x =  w * (T * x + diag(1./diag(M)) * b - previous_x) + previous_x;

end

