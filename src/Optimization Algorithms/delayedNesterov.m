function [ x , parameters] = delayedNesterov( x , previous_x , grad, parameters)
%DELAYEDNESTEROV Applies a novel Nesterov-like method for a supplied alpha
% and beta parameters following the idea of delaying given in [1].
%   Executes y(k+1) = x(k) - alpha grad(x);
%            x(k+1) = (1 + beta) y(k+1) - beta x(k-1)

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

% Nesterov iteration
x = x - parameters.alpha * grad(x);
x = (1 + parameters.beta) * x - parameters.beta * previous_x;

end

