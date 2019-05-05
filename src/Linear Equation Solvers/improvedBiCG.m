function [ x , parameters ] = improvedBiCG( x , ~ , A , b, parameters, ~)
%improvedBiCG Solves the linear equation A x = b applying the Improved Biconjugate Gradient method 
%   for any matrix A. This implementation is highly parallelizable.
% For more information see [11].

% If in the first iteration
if isempty(parameters)
    % Initialize all variables
    parameters.k = 1;
    
    parameters.p = 0;
    parameters.q = 0;
    
    parameters.tilde_v = b - A * x;
    parameters.tilde_w = b - A * x;
    
    parameters.gamma = zeros(1,2);
    parameters.xi = 0;
    
    parameters.tau = 1;
    parameters.rho = ones(2,1);
    
    parameters.kappa = -1;
    
    parameters.a = A' * parameters.q;
    parameters.b = A * parameters.p;
    
    parameters.s = A' * parameters.tilde_w;
    
    parameters.gamma(2) = norm(parameters.tilde_v);
    parameters.xi(2) = norm(parameters.tilde_w);
    
    parameters.rho(2) = parameters.tilde_w' * parameters.tilde_v;
    parameters.epsilon = parameters.s' * parameters.tilde_v;
end

% Load variables
previous_gamma = parameters.gamma(1);
gamma = parameters.gamma(2);

previous_xi = parameters.xi(1);
xi = parameters.xi(2);

previous_rho = parameters.rho(1);
rho = parameters.rho(2);

previous_tau = parameters.tau;

epsilon = parameters.epsilon;

kappa = parameters.kappa;

p = parameters.p;
q = parameters.q;

tilde_v = parameters.tilde_v;
tilde_w = parameters.tilde_w;

s = parameters.s;

% Algorithm
mu = (previous_gamma * previous_xi * rho)/(gamma * previous_tau * previous_rho);

tau = epsilon/rho - gamma * mu;

kappa = -gamma/tau * kappa;

p = 1/gamma * tilde_v - mu *p;

q = 1/xi * s - (gamma * mu)/xi * q;

a = A' * q;

b = A * p;

s = a - tau / xi * s;

tilde_v = b - tau / gamma * tilde_v;

tilde_w = q - tau / xi * tilde_w;

next_gamma = tilde_v' * tilde_v;

next_xi = tilde_w' * tilde_w;

next_rho = tilde_w' * tilde_v;

epsilon = s' * tilde_v;

x = x + kappa * p;

% Store variables
parameters.gamma = [gamma next_gamma];
parameters.xi = [xi next_xi];
parameters.rho = [rho next_rho];

parameters.tau = tau;
parameters.epsilon = epsilon;

parameters.kappa = kappa;

parameters.p = p;
parameters.q = q;

parameters.s = s;
parameters.tilde_v = tilde_v;
parameters.tilde_w = tilde_w;

