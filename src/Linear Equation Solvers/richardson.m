function [ x , parameters] = richardson( x , ~ , A , b, parameters, ~)
%RICHARDSON Solves the linear equation A x = b applying the richardson method 
%   Executes x(k+1) = x(k) + w ( b - A x(k) ) 
%   for w a step parameter.

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

% Richardson iteration
x =  x + parameters.w * ( b - A * x );
end

