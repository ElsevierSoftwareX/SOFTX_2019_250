function w = optWJ(m,L)
%OPTWJ Returns the optimal w that guarantees the fastest convergence for
%the Weighted Jacobi method. 
%   Solves the optimization problem:
%       min    max{ |(1-w)-w*m| , |(1-w)-w*L| }
%   which is equivalent to minimizing the spectral radius of matrix 
%   T_WJ = -w R + (1-w) I, assuming A = I + R.
%   Given the format of the problem, the solution is given by:
%   wopt = 2/(m+L), for m = min(eig(A)) and L = max(eig(A)).

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

% ==== Could solve this optimization ====
% w = sdpvar;
% 
% f1 = abs((1-w)-w*m);
% f2 = abs((1-w)-w*L);
% 
% F = [f1 <= 1,f2 <= 1];
% 
% optimize(F,max(f1,f2));
% 
% w = double(w);

% ==== Closed-form solution ====
w = 2/(m + L);

end

