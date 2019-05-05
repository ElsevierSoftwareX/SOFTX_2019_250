function [ x , parameters] = coordinateDescent( x , ~ , A , b, parameters, ~)
%COORDINATEDESCENT Solves the linear equation A x = b applying the Coordinate Descent. 
%   Executes x(k+1) = x(k) + ( b_i - a_i' * x(k) )/a_ii) * e_i
%   cycling through all rows i. 
%   See [6] [Strohmer] A Randomized Kaczmarz Algorithm with Exponential.

if ~isfield(parameters,'i')
    % Initialize which row to project first
    parameters.i = 0;
end

% % Load current row number
% i = parameters.i;
% 
% % Update to the next row and store it to the next iteration
% i = mod(i , size(A,1)) + 1;
% parameters.i = i;

for i = 1:length(A)
    % Coordinate Descent iteration
    x = x + (b(i,1) - A(i,:) * x) / A(i,i) * e(i,length(x));
end
end

