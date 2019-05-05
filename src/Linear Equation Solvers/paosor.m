function [ x , parameters] = paosor( x , ~ , A , b, parameters, ~)
%PAOSOR Solves the linear equation A x = b applying the Practical Asymptotical Optimal SOR method 
%   Executes x(k+1) = x(k) + wk (I - wk L)^-1 rk
%   for an optimized wk parameter and the residual rk. For more information
%   Check the paper [2] [Meng] A practical asymptotical optimal SOR method

if ~isfield(parameters, 'L')
    % Normalize the system such that A diagonal becomes the I matrix
    b = b ./ diag(A);
    A = A ./ diag(A);
    
    % Store the modified values to compute further residuals
    parameters.A = A;
    parameters.b = b;
    
    % Get the partitioning of A
    L = -tril(A,-1);
    
    % Store matrix L;
    parameters.L = L;
    
    % Initialize rk and wk
    parameters.rk = 1e5;
    parameters.wk = 1;
else
    % Initialize L, A and b with the values in the parameters
    L = parameters.L;
    A = parameters.A;
    b = parameters.b;
end

% Get the previous residual and parameter
previous_wk = parameters.wk;

% Compute residual rk
rk = b - A * x;

% If matrix is positive definite 
[~ , p] = chol(A);
if p == 0 && norm(A-A') < eps
    alpha0 = rk' * rk;
    alpha1 = 2 * rk' * L * rk - rk' * A * rk;
    alpha2 = 3 * rk' * L^2 * rk - 3 * rk' * A * L * rk;
    alpha3 = 4 * rk' * L^3 * rk - 4 * rk' * A * L^2 * rk - 2 * rk' * L' * A * L * rk;
    
    [ w, ~ ] = fullnewton( @(w) 1 + alpha1/alpha0 * w + alpha2/alpha0 * w^2 + alpha3/alpha0 * w^3', @(w) alpha1/alpha0 + alpha2/alpha0 * 2 * w + alpha3/alpha0 * 3 * w^2, 0, 1e-8, 10000 );
    wk = w(end);
    
else
    % wk for a non-symmetric matrix
    beta0 = rk' * A * rk;
    beta1 = 2 * rk' * A * L *rk - rk' * A' * A * rk;
    beta2 = 3 * (rk' * A * L^2 * rk - rk' * A' * A * L * rk);
    beta3 = 4 * rk' * A * L^3 * rk - 4 * rk' * A' * A * L^2 *rk - 2 * rk' * L' * A' * A * L *rk;
    beta4 = 5 * (rk' * A * L^4 * rk - rk' * A' * A * L^3 * rk - rk' * L' * A' * A * L^2 * rk);
    
    [ w, ~ ] = fullnewton(@(w) 1 + beta1/beta0 * w + beta2/beta0 * w^2 + beta3/beta0 * w^3 + beta4/beta0 * w^4, @(w) beta1/beta0 + beta2/beta0 * 2 * w + beta3/beta0 * 3 * w^2 + beta4/beta0 * 4 * w^3, 0, 1e-8, 10000 );
    wk = w(end);
    
    if wk > 2 || wk < 0 || 1 + beta1/beta0 * wk + beta2/beta0 * wk^2 + beta3/beta0 * wk^3 + beta4/beta0 * wk^4 > 1e-6
        wk = previous_wk;
    end
    
end
% end

% PAOSOR iteration
x = x + wk * (eye(length(x)) - wk * L) \ rk;

% Store the current residual and wk for the next iteration
parameters.wk = wk;
