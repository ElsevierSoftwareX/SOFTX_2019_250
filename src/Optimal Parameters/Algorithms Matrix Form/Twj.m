function T = Twj(A,~,w)
%TWJ Computes the transition matrix Twj for the Weighted Jacobi algorithm.
%   Computes :
%   T = I-wD^-1 A

D = diag(diag(A));

T = eye(size(A)) - w * diag(1./diag(D)) * A;

end

