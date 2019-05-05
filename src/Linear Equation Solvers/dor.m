function [ x , parameters] = dor( x , previous_x , A , b, parameters, ~)
%DOR Solves the linear equation A x = b applying the Delayed
%Over-relaxation algorithm to a traditional method such as Jacobi,
%Richardson, Gauss-Seidel, etc. For that it should be supplied the G matrix
%and vector f such that x(k+1) = G x(k) + f.
%   Executes y(k+1) = G x(k) + f
%            x(k+1) = w y(k+1) + (1-w) x(k-1)
%   for a relaxation parameter w. For more details see [1] [Antuono] Delayed Over-Relaxation for Iterative Methods.

if isfield(parameters, 'method') && isfield(parameters, 'parameters')
    % Get the iteration from an already defined method
    [y, ~] = parameters.method(x, 0 , A , b , parameters.parameters, 0);
else
    % Matrix G and vector f of the general function of the algorithm are
    % used to apply the traditional method
    y = parameters.G * x + parameters.f;
end

% Apply the Delayed Over-Relaxation
x = parameters.w * y + (1 - parameters.w) * previous_x;
