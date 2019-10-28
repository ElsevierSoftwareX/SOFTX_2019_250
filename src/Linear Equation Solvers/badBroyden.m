function [ x , parameters] = badBroyden( x , ~ , A , b, parameters, ~)
%BROYDEN Solves the linear equation A x = b applying the "bad" Broyden method
%which is a quasi-Newton method
%   Executes x(k+1) = x(k) - J_k^-1 (A x(k) - b)
%            J_k^-1 = J_(k-1)^-1 + (\Delta x_k - J_(k-1)^-1 \Delta
%            f_k)/(|\Delta f_k |^2) * \Delta f_k'
% 
%   where \Delta x_k = x(k) - x(k-1) and \Delta f_k = A * \Delta x_k.

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

% Broyden iteration
if isempty(parameters)
    % Initialize invJk
    parameters.invJk = eye(length(x));
end

invJk = parameters.invJk;

% Compute step and function evaluation
delta_x = - invJk * (A * x - b);
delta_f = A * delta_x;

% Update estimate for x
x = x + delta_x;

% Update invJk
parameters.invJk = invJk + (delta_x - invJk * delta_f)/(delta_f' * delta_f)* delta_f';
