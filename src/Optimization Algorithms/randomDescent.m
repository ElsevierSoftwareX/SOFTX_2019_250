function [ x , parameters] = randomDescent( x , ~ , grad, parameters)
%RANDOMDESCENT Applies the gradient descent method for a random alpha
%parameter from uniform distribution.
%   Executes x(k+1) = x(k) - alpha grad(x)

% Gradient descent iteration
x = x - 0.8*rand * grad(x);

end

