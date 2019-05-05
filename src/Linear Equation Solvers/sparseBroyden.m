function [ x , parameters] = sparseBroyden( x , ~ , A , b, parameters , max_iterations)
%SPARSEBROYDEN Solves the linear equation A x = b applying the sparse Broyden method
%which is a quasi-Newton method
%   Executes x(k+1) = x(k) - J_k^-1 (A x(k) - b)
%            J^-1_(k+1) = \Pi_(j=0)^(k-1) I + s_(j+1)s_j'/||s_j||^2
%
% The implementation follows the book Iterative Methods for Linear and
% Nonlinear Equations by C. T. Kelley

% Allocate the parameter structure and initialize variable s for the
% direction of improvement
if isempty(parameters)
    parameters.s = zeros(length(x),max_iterations);
    parameters.n = 0;
    
    parameters.s(:,1) = b - A * x;
end

% Use local variables for readibility

s = parameters.s;
n = parameters.n;
    
% Increment iteration counting. It is required for computing correctly the inverse
n = n + 1;

% Take a step towards a better solution
x = x + s(:,n);

% Compute residual
z = b - A * x;

% In this algorithm given the division it is essential to ensure the
% precision of Matlab is not exceeded

if norm(z) < eps
    parameters.s = s;
    parameters.n = n;
    return;
end

% Computing the sparse inverse
for j = 1:(n-1)
    z = z + s(:,j+1)*s(:,j)'*z/(s(:,j)' * s(:,j));
end

% Update the next direction for improvement
s(:,n+1) = z/(1 - s(:,n)'*z/(s(:,n)' * s(:,n)));

% Update the parameters variables
parameters.s = s;
parameters.n = n;