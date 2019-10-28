function [ x , parameters] = kaczmarz( x , ~ , A , b, parameters, ~)
%KACZMARZ Solves the linear equation A x = b applying the Kaczmarz method. 
%   Executes x(k+1) = x(k) + ( b_i - a_i' * x(k) )/(a_i' * a_i) * a_i
%   cycling through all rows i. 
%   See [6] [Strohmer] A Randomized Kaczmarz Algorithm with Exponential.

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

if ~isfield(parameters,'i')
    % Initialize which row to project first
    parameters.i = 0;
end

% % Load current row number
% i = parameters.i;
% 
% % Update to the next row and store it to the next iteration
% i = mod(i , size(A,1)) + 1;
% parameters.i = i;

for i = 1:length(A)
    % Kaczmarz iteration
    x = x + (b(i,1) - A(i,:) * x) / (A(i,:) * A(i,:)') * A(i,:)';
end
end

