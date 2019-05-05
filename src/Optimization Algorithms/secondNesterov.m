function [ x , parameters] = secondNesterov( x , ~ , grad, parameters)
%SECONDNESTEROV Applies the Second Nesterov method for a supplied alpha parameter.
%   Executes v = x(k-1) + (k-2)/(k+1) (x(k-1) - x(k-2)) 
%            x(k) = v - alpha grad(v) 
% See [10] https://www.robots.ox.ac.uk/~vgg/rg/slides/fgrad.pdf
% Equivalent to FISTA but worse behavior if we add projections after each
% update.

if ~isfield(parameters, 'k')
    parameters.k = 1;
    parameters.y = x;
    parameters.v = x;
end

% theta_k
theta_k = 2/(parameters.k+1);

% Gradient descent
parameters.v = parameters.v - parameters.alpha / theta_k * grad( parameters.y );

x = (1-theta_k) * x + theta_k * parameters.v;

parameters.k = parameters.k + 1;
theta_k = 2/(parameters.k+1);

parameters.y = (1 - theta_k) * x + theta_k * parameters.v;
