function [ x , parameters] = broyden( x , ~ , A , b, parameters, ~)
%BROYDEN Solves the linear equation A x = b applying the Broyden method
%which is a quasi-Newton method
%   Executes x(k+1) = x(k) - J_k^-1 (A x(k) - b)
%            J_k^-1 = J_(k-1)^-1 + (\Delta x_k - J_(k-1)^-1 \Delta f_k)/
%            (\Delta x_k' J_(k-1)^-1 \Delta f_k) \Delta x_k' J_(k-1)^-1
% 
%   where \Delta x_k = x(k) - x(k-1) and \Delta f_k = A * \Delta x_k.

% Broyden iteration
if isempty(parameters)
    % Initialize invJk
    parameters.invJk = eye(length(x));
end

invJk = parameters.invJk;

% Compute step and function evaluation
delta_x = - invJk * (A * x - b);
delta_f = A * delta_x;

% Update estimate for x
x = x + delta_x;

% Update invJk
parameters.invJk = invJk + (delta_x - invJk * delta_f)/(delta_x' * invJk * delta_f)* delta_x' * invJk;
