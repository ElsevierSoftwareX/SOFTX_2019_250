function [ x , parameters] = momentum( x , ~ , grad, parameters)
%MOMENTUM Applies the Momentum method for a supplied alpha
% and beta parameters
%   Executes x(k+1) = x(k) - alpha grad(x(k)) + beta v(k)
%            v(k+1) = - alpha grad(x(k)) + beta v(k) 
% See https://jlmelville.github.io/mize/nesterov.html. 
% If we replace the update of v(k+1) for:
%            v(k+1) = - alpha grad(x(k+1)) + beta v(k)
% we obtain a smoother algorithm equal to nesterov's without any
% projections but worse when those are present.

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

if ~isfield(parameters, 'v')
    parameters.v = zeros(length(x),1);
end

% Compute step 
d = - parameters.alpha * grad(x) + parameters.beta * parameters.v;

% Momentum iteration
x = x + d;
% Velocity
parameters.v = d;

end

