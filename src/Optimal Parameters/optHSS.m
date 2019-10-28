function [alpha] = optHSS(A, method)
%OPTHSS Estimates the optimal alpha parameter for HSS iteration.
%   Sets alpha to be the root of:
%   4n alpha^3 + 3 eta alpha^2 + 2 w alpha + gamma with:
%       eta = -2 tr(H), w = tr(H^2) - tr(S^2) and gamma = 2 tr(HS^2).
%   Method 1 - See [5] [Huang] A practical formula for computing optimal parameters
%   in the HSS iteration methods
%   Method 2 - See [7] [Chen] On choices of iteration parameter in HSS method

% Copyright 2019 Daniel Silvestre
% This file is part of OPTool.
%
% OPTool is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% OPTool is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with OPTool.  If not, see <https://www.gnu.org/licenses/>.

% Build H and S
H = 0.5 * (A + A');
S = 0.5 * (A - A');

switch method 
    case 1 % [5]
        % Define the scalars of the polynomial to solve
        eta = -2 * trace(H);
        w = trace(H^2) - trace(S^2);
        gamma = 2 * trace(H*S^2);
        
        % Solve using the Newton method the polynomial
        [ alpha, ~ ] = fullnewton( @(x) 4*length(A)*x^3 + 3*eta*x^2 + 2*w*x + gamma,@(x) 12*length(A)*x^2 + 6*eta*x + 2*w, 1e2, 1e-8, 1e5 );
        alpha = alpha(end);
    case 2
        lmin = min(eig(H));
        lmax = max(eig(H));
        
        smax = max(svd(S));
        smin = min(svd(S));
        
        % Solve using the Newton method the polynomial
        [ alpha, ~ ] = fullnewton( @(x) 2*(lmax-lmin)*x^3 + (lmax^2-lmin^2-(smax^2-smin^2))*x^2 + 2*(smin^2*lmax-smax^2*lmin)*x - smax^2*lmin^2,@(x) 6*(lmax-lmin)*x^2 + 2*(lmax^2-lmin^2-(smax^2-smin^2))*x + 2*(smin^2*lmax-smax^2*lmin), 1e2, 1e-8, 1e5 );
        alpha = alpha(end);
        
    otherwise
        alpha = 1;
end

