function [ x , parameters] = richardson( x , ~ , A , b, parameters, ~)
%RICHARDSON Solves the linear equation A x = b applying the richardson method 
%   Executes x(k+1) = x(k) + w ( b - A x(k) ) 
%   for w a step parameter.

% Richardson iteration
x =  x + parameters.w * ( b - A * x );
end

