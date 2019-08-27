function [ x , parameters] = nesterovLTV( x , previous_x , grad, parameters)
%NESTEROVLTV Applies the Nesterov method for a supplied alpha parameter and
%using beta = (k-1)/(k+2)
%   Executes x(k+1) = y(k) - alpha grad(y);
%            y(k) = (1 + beta) x(k) - beta x(k-1)

% Access time variable k to compute beta
if ~isfield(parameters,'k')
    k = 0;
    symx = sym('x',[length(x),1]);
    parameters.alpha = 1/max(eig(double(jacobian(grad(symx),symx))));
else
    k = parameters.k + 1;
end

% Compute parameter beta
parameters.beta = (k-1)/(k+2);

% Compute auxiliary variable y(k)
y = (1 + parameters.beta) * x - parameters.beta * previous_x;

% Nesterov iteration
x = y - parameters.alpha * grad(y);

% Store time variable k
parameters.k = k;
end

