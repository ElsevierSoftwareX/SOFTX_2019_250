function w = optWJ(m,L)
%OPTWJ Returns the optimal w that guarantees the fastest convergence for
%the Weighted Jacobi method. 
%   Solves the optimization problem:
%       min    max{ |(1-w)-w*m| , |(1-w)-w*L| }
%   which is equivalent to minimizing the spectral radius of matrix 
%   T_WJ = -w R + (1-w) I, assuming A = I + R.
%   Given the format of the problem, the solution is given by:
%   wopt = 2/(m+L), for m = min(eig(A)) and L = max(eig(A)).

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

