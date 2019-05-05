function [ x , parameters] = mrDOR( x , previous_x , A , b, parameters, ~)
%MRDOR Solves the linear equation A x = b applying the Delayed
%Over-relaxation algorithm to a traditional Richardson method and using the
%idea of minimal residual to compute the values for w and t.
%   Executes y(k+1) = x(k) + t * (A x - b)
%            x(k+1) = w y(k+1) + (1-w) x(k-1)
%   for a relaxation parameter w and step t. For more details see [1] [Antuono] Delayed Over-Relaxation for Iterative Methods.

% Compute tk for Richardson iteration based on the residual of x at time k,
% rk
rk = A * x - b;
% Product A * rk to save computation
A_rk = A * rk;
% Paramter tk definition
tk = -(rk' * A_rk) / (A_rk' * A_rk);

% Richardson iteration
y =  x + tk * rk;

% Compute Over-Relaxation parameter w from r_(k-1) and residual for y
residual_previous_k = A * previous_x - b;
residual_y = A * y - b;
% Subtraction between r_(k-1) - y(k + 1) to save computation
residual_rk1_y = residual_previous_k - residual_y;
% Parameter wk definition
wk = (residual_previous_k' * residual_rk1_y ) / ( residual_rk1_y' * residual_rk1_y);

% Apply the Delayed Over-Relaxation
x = wk * y + (1 - wk) * previous_x;
