function [ x , parameters] = jacobi( x , ~ , A , b, parameters, ~)
%JACOBI Solves the linear equation A x = b applying the jacobi method 
%   Executes x(k+1) = D^-1 ( b - R x(k) ) 
%   for A = D + R with D a diagonal matrix.

% Partition matrix A = D + R
D = diag(diag(A));
R = A - D;

% Residual vector
residual = b - R * x;

% Jacobi iteration
x =  diag(1./diag(D)) * residual;

end

