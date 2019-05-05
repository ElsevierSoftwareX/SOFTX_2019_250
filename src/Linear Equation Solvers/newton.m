function [ x , parameters] = newton( ~ , ~ , A , b, parameters, ~)
%NEWTON Solves the linear equation A x = b applying the Newton-Raphson method 
%   Executes x(k+1) = x(k) - A^-1 (A x(k) - b) 

% Newton-Raphson iteration
x =  A\ b;

end

