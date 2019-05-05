function w_max = optRichardson(m,L)
%OPTRICHARDSON Returns the optimal w that guarantees the fastest convergence for
%the Richardson method. 
%   Solves the optimization problem:
%       min    max{ |1-w*m| , |1-w*L| }
%   which is equivalent to minimizing the spectral radius of matrix 
%   T_Richardson = I - w A.
%   Given the format of the problem, the solution is given by:
%   wopt = 2/(m+L).

% ==== Closed-form solution ====
w_max = 2/(m + L);
end

