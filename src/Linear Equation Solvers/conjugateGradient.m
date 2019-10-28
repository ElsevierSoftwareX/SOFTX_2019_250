function [ x , parameters ] = conjugateGradient( x , ~ , A , b, parameters, ~)
%conjugateGradient Solves the linear equation A x = b applying the Conjugate Gradient method 
%   for a positive definite matrix A. 
%   Executes x(k+1) = x(k) + alpha p_k for a new conjugate direction p_k.

% % Conjugate Gradient only applies to symmetric matrices
% if norm(A-A') > 1E-6
%     % Solving the normal equation
%     fprintf('A is not symmetric. Solving the normal equation A\''A x = A\''b\n');
%     A = A'*A;
%     b = A'*b;
% end

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

% If in the first iteration
if isempty(parameters)
    % Allocate the parameters variable
    parameters = zeros(length(x),2);
    
    % Compute r_0
    parameters(:,2) = b - A * x;
    % p_0 = r_0
    parameters(:,1) = parameters(:,2);
end
   
% Compute products used twice. A *direction p_k and r_k' * r_k
Ap = A * parameters(:,1);
normSquaredRk = parameters(:,2)' * parameters(:,2);

% Compute new alpha
alpha = normSquaredRk / (parameters(:,1)' * Ap);

% Iterate using the current direction
x = x + alpha * parameters(:,1);

% Update the residual r_k+1
parameters(:,2) = parameters(:,2) - alpha * Ap;

% Update the new direction p_k+1
parameters(:,1) = parameters(:,2) + ((parameters(:,2)' * parameters(:,2)) / normSquaredRk) * parameters(:,1);

end

