function [ x , parameters] = cbb( x , ~ , grad, parameters)
%CBB Applies the Cauchy-Barzilai-Borwein method for quadratic problems.
%   Executes x(k+1) = x(k) - alpha grad(x) + alpha^2 Q grad(x)
%   where grad(x) = Qx - b and alpha = ||grad(x)||^2 / (grad(x)' Q grad(x))
%   For additional information see [12].

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


% Check if supplied with the Q matrix defining the quadratic problem
if ~isfield(parameters,'Q')
    disp('It is missing the Q matrix defining the quadratic problem. I will assume the identity matrix.');
    parameters.Q = eye(length(x));
end

% Load Q matrix
Q = parameters.Q;

% Compute gradient
gk = grad(x);

% Compute search direction
hk = Q * gk;

% Compute step-size
alpha = (gk' * gk) / (gk' * hk);

% Gradient descent iteration
x = x - 2*alpha *gk +alpha^2 * hk;

end

