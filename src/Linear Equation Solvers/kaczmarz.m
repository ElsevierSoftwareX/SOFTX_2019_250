function [ x , parameters] = kaczmarz( x , ~ , A , b, parameters, ~)
%KACZMARZ Solves the linear equation A x = b applying the Kaczmarz method. 
%   Executes x(k+1) = x(k) + ( b_i - a_i' * x(k) )/(a_i' * a_i) * a_i
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
    % Kaczmarz iteration
    x = x + (b(i,1) - A(i,:) * x) / (A(i,:) * A(i,:)') * A(i,:)';
end
end
