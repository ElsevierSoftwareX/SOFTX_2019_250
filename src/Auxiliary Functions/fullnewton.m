function [ x, ex ] = fullnewton( f, df, x0, tol, nmax )
%
% NEWTON Newton's Method
%   Newton's method for finding successively better approximations to the 
%   zeroes of a real-valued function.
%
% Input:
%   f - input funtion
%   df - derived input function
%   x0 - inicial aproximation
%   tol - tolerance
%   nmax - maximum number of iterations
%
% Output:
%   x - aproximation to root
%   ex - error estimate
%
% Example:
%	[ x, ex ] = newton( 'exp(x)+x', 'exp(x)+1', 0, 0.5*10^-5, 10 )
%
% Author:	Tashi Ravach
% Revised:  Daniel Silvestre
% Version:	1.0
% Date:     16/04/2007
%

    if nargin == 3
        tol = 1e-4;
        nmax = 1e1;
    elseif nargin == 4
        nmax = 1e1;
    elseif nargin ~= 5
        error('newton: invalid input parameters');
    end
    
    x = zeros(1,nmax);
    ex = zeros(1,nmax);
    
    x(1) = x0 - (f(x0)/df(x0));
    ex(1) = norm(f(x(1)));
    k = 2;
    while (ex(k-1) >= tol) && (k <= nmax)
        x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
        ex(k) = norm(f(x(k)));
        k = k+1;
    end
    x = x(1:k-1);
end