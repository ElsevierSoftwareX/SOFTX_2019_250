function [ x , parameters ] = cgne( x , ~ , A , b, parameters, ~)
%CGNE Solves the linear equation A x = b applying the Conjugate Gradient to the
% Normal equation.
%   Computes directions conjugate in which to reduce the residual. 
%   See pg 241 of [8] [Kyrchei] Advances in Linear Algebra Research.

% If in the first iteration
if isempty(parameters)
    % Compute the initial residual
    parameters.r = b - A * x;
    
    % Compute initial direction 
    parameters.p = A' * parameters.r;
end

% Load variables
r = parameters.r;
p = parameters.p;

% Compute step size for the direction p
alpha = (r' * r) / (p' * p);

% Update x with direction p and step size alpha
x = x + alpha *p;

% Compute next residual
next_r = r - alpha * A * p;

% Step size to obtain the new direction
beta = (next_r' * next_r) / (r' * r);

% Get new direction
p = A' * next_r + beta * p;

% Save residual and direction for next iteration
parameters.p = p;
parameters.r = next_r;