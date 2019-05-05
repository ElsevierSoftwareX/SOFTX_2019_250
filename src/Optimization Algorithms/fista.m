function [ x , parameters] = fista( x , previous_x , grad, parameters)
%FISTA Applies the Fast Iterative Soft-thresholding Algorithm (FISTA)
%method for a supplied alpha parameter.
%   Executes v = x(k-1) + (k-2)/(k+1) (x(k-1) - x(k-2)) 
%            x(k) = v - alpha grad(v) 
% See [9] https://www.cs.cmu.edu/~ggordon/10725-F12/slides/09-acceleration.pdf

if ~isfield(parameters, 'k')
    parameters.k = 1;
end

% Velocity
v = x + (parameters.k-2)/(parameters.k+1) * (x - previous_x);

% Gradient Descent
x = v - parameters.alpha * grad( v );

% Update k
parameters.k = parameters.k + 1;
end

