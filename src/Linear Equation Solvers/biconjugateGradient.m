function [ x , parameters ] = biconjugateGradient( x , ~ , A , b, parameters, ~)
%conjugateGradient Solves the linear equation A x = b applying the Biconjugate Gradient method 
%   for any matrix A. 

% If in the first iteration
if isempty(parameters)
    % Allocate the parameters variable
    parameters = zeros(length(x),2);
    
    % Compute r_0
    parameters(:,1) = b - A * x;
    % Compute hat_r_0
    parameters(:,2) = parameters(:,1);
    % Compute rho_0, alpha, w0 = 1
    parameters(1:3,3) = ones(3,1);
    % Compute v_0, p0 = 0
    parameters(:,4:5) = zeros(length(x),2);
end

% Compute rho_i and store rho_{i-1}
old_rho = parameters(1,3);
parameters(1,3) = parameters(:,2)' * parameters(:,1);

% Compute beta
beta = (parameters(1,3)/old_rho) * (parameters(2,3)/parameters(3,3));

% Compute new p_i direction
parameters(:,5) = parameters(:,1) + beta * (parameters(:,5) - parameters(3,3) * parameters(:,4));

% Compute new v_i
parameters(:,4) = A * parameters(:,5);

% Compute alpha
parameters(2,3) = parameters(1,3) / (parameters(:,2)' * parameters(:,4));

% Compute partial new value of x 
h = x + parameters(2,3) * parameters(:,5);

% Correction term
s = parameters(:,1) - parameters(2,3) * parameters(:,4);
t = A * s;

% New w_k
parameters(3,3) = (t' * s) / (t' * t);

% New value of x
x = h + parameters(3,3) * s;

% Next direction for future iteration
parameters(:,1) = s - parameters(3,3) * t;

end

