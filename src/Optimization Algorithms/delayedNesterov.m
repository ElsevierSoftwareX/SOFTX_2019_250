function [ x , parameters] = delayedNesterov( x , previous_x , grad, parameters)
%DELAYEDNESTEROV Applies a novel Nesterov-like method for a supplied alpha
% and beta parameters following the idea of delaying given in [1].
%   Executes y(k+1) = x(k) - alpha grad(x);
%            x(k+1) = (1 + beta) y(k+1) - beta x(k-1)

% Nesterov iteration
x = x - parameters.alpha * grad(x);
x = (1 + parameters.beta) * x - parameters.beta * previous_x;

end

