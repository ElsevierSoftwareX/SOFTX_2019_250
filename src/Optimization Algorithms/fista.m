function [ x , parameters] = fista( x , previous_x , grad, parameters)
%FISTA Applies the Fast Iterative Soft-thresholding Algorithm (FISTA)
%method for a supplied alpha parameter.
%   Executes v = x(k-1) + (k-2)/(k+1) (x(k-1) - x(k-2)) 
%            x(k) = v - alpha grad(v) 
% See [9] https://www.cs.cmu.edu/~ggordon/10725-F12/slides/09-acceleration.pdf

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

if ~isfield(parameters, 'k')
    parameters.k = 1;
end

% Velocity
v = x + (parameters.k-2)/(parameters.k+1) * (x - previous_x);

% Gradient Descent
x = v - parameters.alpha * grad( v );

% Update k
parameters.k = parameters.k + 1;
end

