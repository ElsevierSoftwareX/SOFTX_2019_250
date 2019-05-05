function [ x , parameters] = weightedJacobi( x , ~ , A , b, parameters, ~)
%WEIGHTEDJACOBI Solves the linear equation A x = b applying the relaxed jacobi method 
%   Executes x(k+1) = w D^-1 ( b - R x(k) ) + (1-w) x(k)
%   for A = D + R with D a diagonal matrix.

% Partition matrix A = D + R
D = diag(diag(A));
R = A - D;

% Residual vector
residual = b - R * x;

% Jacobi iteration
x =  parameters.w * diag(1./diag(D)) * residual + (1-parameters.w) * x;

end

