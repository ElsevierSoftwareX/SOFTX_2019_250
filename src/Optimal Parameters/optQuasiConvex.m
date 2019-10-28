function w_max = optQuasiConvex(A,functionToGetTransition,G)
%OPTQUASICONVEX Returns the optimal w that guarantees the fastest convergence for
%any algorithm with quasiconvex property on the optimization of T that might depend on G.
%   Solves the optimization problem:
%       min    max{ lambda( T ) }
%   Given the format of the problem, the solution is found through
%   bissetion.

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

lower = 0;
upper = 2;
epsilon = 1E-5;

while (upper-lower) >= epsilon
    w = (lower+upper)/2;
    
    Tw = functionToGetTransition(A,G,w);
    dTw = functionToGetTransition(A,G,w+epsilon);
    
    if spectralRadius(dTw) < spectralRadius(Tw)
        lower = w+epsilon;
    else
        upper = w;
    end
end
w_max = (lower+upper)/2;
end

