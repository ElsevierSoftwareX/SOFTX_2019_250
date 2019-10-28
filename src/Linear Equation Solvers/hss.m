function [ x , parameters] = hss( x , ~ , A , b, parameters, ~)
%HSS Solves the linear equation A x = b applying the HSS method
%   Executes y(k+1) = (alpha I + H)^-1 *( (alpha I - S) * x(k) + b)
%            x(k+1) = (alpha I + S)^-1 *( (alpha I - H) * y(k+1) + b)
%   See [5] A practical formula for computing optimal parameters in the HSS iteration methods.

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

if ~isfield(parameters, 'H') || ~isfield(parameters, 'S')
    % Partition matrix A = H + S
    H = 0.5 * (A + A');
    S = 0.5 * (A - A');
    
    % Store all variables for the next iteration
    parameters.H = H;
    parameters.S = S;
end

% Load variables
H = parameters.H;
S = parameters.S;
alpha = parameters.alpha;

% HSS iteration
y = (alpha * eye(length(x)) + H) \ ( (alpha * eye(length(x)) - S) * x + b);
x = (alpha * eye(length(x)) + S) \ ( (alpha * eye(length(x)) - H) * y + b);

end

