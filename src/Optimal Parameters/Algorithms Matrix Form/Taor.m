function T = Taor(A,~,w,r)
%TAOR Computes the transition matrix Taor for the AOR algorithm.
%   Computes :
%   T = (I - r * L)^-1 ((1 - w) * I + (w - r) * L + w * U);

% Compute the partition of matrix A
D = diag(diag(A));
A_U = -triu(A,1);
A_L = -tril(A,-1);
% Compute the U and L matrices
L = diag(1./diag(D)) * A_L;
U = diag(1./diag(D)) * A_U;

T = (eye(length(A)) - r * L) \ ((1 - w) * eye(length(A)) + (w - r) * L + w * U);

end

