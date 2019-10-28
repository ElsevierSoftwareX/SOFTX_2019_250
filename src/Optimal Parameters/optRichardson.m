function w_max = optRichardson(m,L)
%OPTRICHARDSON Returns the optimal w that guarantees the fastest convergence for
%the Richardson method. 
%   Solves the optimization problem:
%       min    max{ |1-w*m| , |1-w*L| }
%   which is equivalent to minimizing the spectral radius of matrix 
%   T_Richardson = I - w A.
%   Given the format of the problem, the solution is given by:
%   wopt = 2/(m+L).

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

% ==== Closed-form solution ====
w_max = 2/(m + L);
end

