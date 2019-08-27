function T = Tj(A,~,~)
%TSOR Computes the transition matrix Tj for the Jacobi algorithm.
%   Computes :
%   T = -D^-1*R 

D = diag(diag(A));
R = A-D;

T = -inv(D)*R;

end

