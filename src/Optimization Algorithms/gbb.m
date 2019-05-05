function [ x , parameters] = gbb( x , ~ , grad, parameters)
%CBB Applies the Global Barzilai-Borwein method for quadratic problems.
%   Executes x(k+1) = x(k) - alpha grad(x)
%   where alpha is the result of a nonmonotone search.
%   For additional information see [13].

% First iteration
if ~isfield(parameters,'fk')
    % Check if it is supplied the appropriate parameters
    if ~isfield(parameters,'M')
        disp('M was not specified. I will use the default value of 10.');
        parameters.M = 10;
    end
    if ~isfield(parameters,'rho')
        disp('rho was not specified. I will use the default value of 0.5.');
        parameters.rho = 0.5;
    end
    if ~isfield(parameters,'delta')
        disp('delta was not specified. I will use the default value of 0.5.');
        parameters.delta = 0.5;
    end
    if ~isfield(parameters,'sigma')
        disp('sigma was not specified. I will use the default value of 0.8.');
        parameters.sigma = 0.8;
    end
    if ~isfield(parameters,'amin') && ~isfield(parameters,'amax')
        disp('The interval [amin amax] was not specified. I will use the default [0.1 2.0].');
        parameters.amin = 0.1;
        parameters.amax = 2;
    end
    
    % Allocate space for the evaluation of the gradient vector
    if ~isfield(parameters,'f')
        disp('It is missing the function f');
        return;
    else
        parameters.fk = -inf(1,parameters.M);
    end
    
    if ~isfield(parameters,'alpha')
        parameters.alpha = parameters.delta;
    end
    parameters.g = grad(x);
end

% Load parameters
alpha = parameters.alpha;
amin = parameters.amin;
amax = parameters.amax;
fk = parameters.fk;
rho = parameters.rho;
f = parameters.f;
sigma = parameters.sigma;
previous_g = parameters.g;
delta = parameters.delta;

% Test if current alpha is within the acceptable bounds
if alpha <= amin || alpha >= amax
    alpha = delta;
end

lambda = 1/alpha;

% Compute gradient
g = grad(x);

% Update fk
fk = [fk(2:end) f(x)];

while f(x - lambda * g) > max(fk - rho * lambda * (g' * g))
    lambda = sigma * lambda;
end

% Gradient descent iteration
x = x - lambda * g;

% Update alpha
y = g - previous_g;
parameters.alpha = - (g' * y) / (lambda * (g' * g));

% Update fk values 
parameters.fk = fk;
parameters.g = g;