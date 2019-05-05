function [ x , parameters] = gaussSeidel( x , ~ , A , b, parameters, ~)
%gaussSeidel Solves the linear equation A x = b applying the Gauss-Seidel method 
%   Executes x(k+1) = L^-1 ( b - U x(k) ) 
%   for A = L + U with L a lower triangular matrix and U a strict upper
%   triangular matrix. To avoid the inverse, the element-wise definition is
%   used: 
%   x_i(k+1) = 1/a_ii (b_i - \sum_{j=1}^{i-1} a_ij x_j(k+1) - \sum_{j=i+1}^{n} a_ij x_j(k))

n = length(x);

for i = 1:n
    % Gauss-Seidel element-wise iteration
    x(i) = (1/A(i, i))*( b(i) - A(i, [1:i-1 i+1:end])*x([1:i-1 i+1:end],1) );
end

