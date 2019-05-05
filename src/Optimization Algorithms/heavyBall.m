function [ x , parameters] = heavyBall( x , previous_x , grad, parameters)
%HEAVYBALL Applies the Heavy-Ball method for a supplied alpha
% and beta parameters
%   Executes x(k+1) = x(k) - alpha grad(x) + beta ( x(k) - x(k-1) )

% Heavy-Ball iteration
x = x - parameters.alpha * grad(x) + parameters.beta * (x - previous_x);

end

