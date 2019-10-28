function [ x , parameters] = aaj( x , previous_x , A , b, parameters, ~)
%AAJ Solves the linear equation A x = b applying the Alternating
%Anderson-Jacobi.
%   Executes x(k+1) = x(k) + Bk * f(x(k))
%   for f(x) = x(k+1) - x(k) using Jacobi
%   and Bk = wI                                       if (k+1)/p \notin N
%            betaI - (Xk + beta Fk) (Fk' * Fk)^-1 Fk' if (k+1)/p \in N.
%   Xk = [(x(k-m+1) - x(k-m)) ... (x(k) - x(k-1)]
%   Fk = [(f(x(k-m+1)) - f(x(k-m))) ... (f(x(k)) - f(x(k-1))].
%   For more information check:
%   [3] [Pask] Anderson acceleration of the Jacobi iterative method
%   An efficient alternative to Krylov methods for large sparse linear systems

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

if ~isfield(parameters,'Xk') || ~isfield(parameters,'Fk')
    % Declare the variables Xk and Fk in the parameters
    parameters.Xk = zeros(length(x),parameters.m);
    parameters.Fk = zeros(length(x),parameters.m); 
    
    % Keep record of current iteration
    parameters.k = 1;
    
    % Partition matrix A as in the Jacobi method
    parameters.D = diag(diag(A));
    parameters.R = A - parameters.D;
end

% Load the parameters and variables
Xk = parameters.Xk;
Fk = parameters.Fk;
k = parameters.k;
beta = parameters.beta;
w = parameters.w;
p = parameters.p;

% Shift Xk and Fk 1 entry to the left
Xk(:,1:end-1) = Xk(:,2:end);
Fk(:,1:end-1) = Fk(:,2:end);

% Update Xk and Fk
Xk(:,end) = x - previous_x;
Fk(:,end) = -A * Xk(:,end);

% Select Bk
if mod(k+1, p) == 0
    % Acceleration step
    Bk = beta * eye(length(x)) - (Xk + beta * Fk)*pinv(Fk' * Fk) * Fk';
else
    % Weighted Jacobi step
    Bk = w * eye(length(x));
end

% AAJ iteration
x = x + Bk * (b - A * x);

% Update time and save all variables
parameters.Xk = Xk;
parameters.Fk = Fk;
parameters.k = k + 1;
end

