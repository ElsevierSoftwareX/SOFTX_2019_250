function [ x , parameters] = hss( x , ~ , A , b, parameters, ~)
%HSS Solves the linear equation A x = b applying the HSS method
%   Executes y(k+1) = (alpha I + H)^-1 *( (alpha I - S) * x(k) + b)
%            x(k+1) = (alpha I + S)^-1 *( (alpha I - H) * y(k+1) + b)
%   See [4] [Wen] Quasi-Chebyshev accelerated iteration methods based on 
%   optimization for linear systems.

if ~isfield(parameters, 'H') || ~isfield(parameters, 'S')
    % Partition matrix A = H + S
    H = 0.5 * (A + A');
    S = 0.5 * (A - A');
    
    % Store all variables for the next iteration
    parameters.H = H;
    parameters.S = S;
end

% Load variables
H = parameters.H;
S = parameters.S;
alpha = parameters.alpha;

% HSS iteration
y = (alpha * eye(length(x)) + H) \ ( (alpha * eye(length(x)) - S) * x + b);
x = (alpha * eye(length(x)) + S) \ ( (alpha * eye(length(x)) - H) * y + b);

end

