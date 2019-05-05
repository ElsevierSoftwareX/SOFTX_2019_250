function [ x , parameters] = momentum( x , ~ , grad, parameters)
%MOMENTUM Applies the Momentum method for a supplied alpha
% and beta parameters
%   Executes x(k+1) = x(k) - alpha grad(x(k)) + beta v(k)
%            v(k+1) = - alpha grad(x(k)) + beta v(k) 
% See https://jlmelville.github.io/mize/nesterov.html. 
% If we replace the update of v(k+1) for:
%            v(k+1) = - alpha grad(x(k+1)) + beta v(k)
% we obtain a smoother algorithm equal to nesterov's without any
% projections but worse when those are present.

if ~isfield(parameters, 'v')
    parameters.v = zeros(length(x),1);
end

% Compute step 
d = - parameters.alpha * grad(x) + parameters.beta * parameters.v;

% Momentum iteration
x = x + d;
% Velocity
parameters.v = d;

end

