function [ x , parameters] = cbb( x , ~ , grad, parameters)
%CBB Applies the Cauchy-Barzilai-Borwein method for quadratic problems.
%   Executes x(k+1) = x(k) - alpha grad(x) + alpha^2 Q grad(x)
%   where grad(x) = Qx - b and alpha = ||grad(x)||^2 / (grad(x)' Q grad(x))
%   For additional information see [12].

% Check if supplied with the Q matrix defining the quadratic problem
if ~isfield(parameters,'Q')
    disp('It is missing the Q matrix defining the quadratic problem. I will assume the identity matrix.');
    parameters.Q = eye(length(x));
end

% Load Q matrix
Q = parameters.Q;

% Compute gradient
gk = grad(x);

% Compute search direction
hk = Q * gk;

% Compute step-size
alpha = (gk' * gk) / (gk' * hk);

% Gradient descent iteration
x = x - 2*alpha *gk +alpha^2 * hk;

end

