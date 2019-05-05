function [ x , parameters] = sor( x , ~ , A , b, parameters, ~)
%SOR Solves the linear equation A x = b applying the Successive Over-relaxation method 
%   Executes x(k+1) = (D+wL)^-1 ( w b - (w U + (w - 1) D) x(k) ) 
%   for A = D + L + U with L a strict lower triangular matrix, U a strict upper
%   triangular matrix and D a diagonal matrix. To avoid the inverse, the element-wise definition is
%   used: 
%   x_i(k+1) = (1-w) x_i(k) + w/a_ii (b_i - \sum_{j=1}^{i-1} a_ij x_j(k+1) - \sum_{j=i+1}^{n} a_ij x_j(k))

n = length(x);

for i = 1:n
    % Successive Over-relaxation element-wise iteration
    x(i) = (1-parameters.w) * x(i) + (parameters.w/A(i, i))*( b(i) - A(i, [1:i-1 i+1:end])*x([1:i-1 i+1:end],1) );
end

