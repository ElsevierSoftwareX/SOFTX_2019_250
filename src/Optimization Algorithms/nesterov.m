function [ x , parameters] = nesterov( x , previous_x , grad, parameters)
%NESTEROV Applies the Nesterov method for a supplied alpha
% and beta parameters
%   Executes x(k+1) = y(k) - alpha grad(y);
%            y(k) = (1 + beta) x(k) - beta x(k-1)

% Compute auxiliary variable y(k)
y = (1 + parameters.beta) * x - parameters.beta * previous_x;

% Nesterov iteration
x = y - parameters.alpha * grad(y);

end

