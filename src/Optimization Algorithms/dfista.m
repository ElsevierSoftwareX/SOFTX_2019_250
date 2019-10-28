function [ x , parameters] = dfista( x , ~ , grad, parameters)
%DFISTA Applies the Descent Fast Iterative Soft-thresholding Algorithm (DFISTA)
%method for a supplied alpha parameter.
%   Executes u(k) = y(k-1) -alpha grad( y(k-1) ) 
%            x(k) = u(k)    if f(u) <= f(x(k-1))
%                   x(k-1)  else
%            v(k) = x(k-1) + 1/theta_k * ( u(k) - x(k-1) )
%            y(k) = (1 - theta_{k+1}) x(k) + theta_{k+1} v(k)
% See [10] https://www.robots.ox.ac.uk/~vgg/rg/slides/fgrad.pdf

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
    parameters.y = x;
end

% Tentative Gradient descent
u = parameters.y - parameters.alpha * grad( parameters.y );

% If u(k) has worst cost than x stay with x
if parameters.f(u) <= parameters.f(x)
    next_x = u;
else
    next_x = x;
end

% theta_k
theta_k = 2/(parameters.k+1);

% v 
v = x + 1/theta_k * (u - x);

% update k 
parameters.k = parameters.k + 1;
theta_k = 2/(parameters.k+1);

% y
parameters.y = (1 - theta_k) * next_x + theta_k * v;

x = next_x;

