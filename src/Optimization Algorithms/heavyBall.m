function [ x , parameters] = heavyBall( x , previous_x , grad, parameters)
%HEAVYBALL Applies the Heavy-Ball method for a supplied alpha
% and beta parameters
%   Executes x(k+1) = x(k) - alpha grad(x) + beta ( x(k) - x(k-1) )

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

% Heavy-Ball iteration
x = x - parameters.alpha * grad(x) + parameters.beta * (x - previous_x);

end

