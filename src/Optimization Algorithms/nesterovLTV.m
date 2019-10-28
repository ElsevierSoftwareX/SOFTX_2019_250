function [ x , parameters] = nesterovLTV( x , previous_x , grad, parameters)
%NESTEROVLTV Applies the Nesterov method for a supplied alpha parameter and
%using beta = (k-1)/(k+2)
%   Executes x(k+1) = y(k) - alpha grad(y);
%            y(k) = (1 + beta) x(k) - beta x(k-1)

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

% Access time variable k to compute beta
if ~isfield(parameters,'k')
    k = 0;
    symx = sym('x',[length(x),1]);
    parameters.alpha = 1/max(eig(double(jacobian(grad(symx),symx))));
else
    k = parameters.k + 1;
end

% Compute parameter beta
parameters.beta = (k-1)/(k+2);

% Compute auxiliary variable y(k)
y = (1 + parameters.beta) * x - parameters.beta * previous_x;

% Nesterov iteration
x = y - parameters.alpha * grad(y);

% Store time variable k
parameters.k = k;
end

