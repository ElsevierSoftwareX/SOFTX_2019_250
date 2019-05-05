function [ x , parameters] = gradientDescent( x , ~ , grad, parameters)
%GRADIENTDESCENT Applies the gradient descent method for a supplied alpha
%parameter
%   Executes x(k+1) = x(k) - alpha grad(x)

% Gradient descent iteration
x = x - parameters.alpha * grad(x);

end

