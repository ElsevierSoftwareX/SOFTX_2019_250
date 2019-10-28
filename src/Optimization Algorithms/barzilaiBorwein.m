function [ x , parameters] = barzilaiBorwein( x , ~ , grad, parameters)
%BARZILAIBORWEIN Applies the Barzilai-Borwein gradient descent method.
%   Executes x(k+1) = x(k) - alpha_k grad(x)
%            alpha_k = s_{k-1}' s_{k-1} / s_{k-1}' y_{k-1}
%   For s_{k-1} = x(k) - x(k-1) and y_{k-1} = grad(x(k)) - grad(x(k-1)).

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


% First iteration
if ~isfield(parameters,'yk') || ~isfield(parameters,'dk') || ~isfield(parameters,'gk')
    % Compute gradient of current estimate
    grad_current_x = grad(x);
    
    % Compute step to update the current estimate
    dk = - 0.5 * grad_current_x;
    
    % Gradient descent iteration
    x = x + dk;
    
    % Store parameters for next update
    parameters.gk = grad(x);
    parameters.yk = parameters.gk - grad_current_x;
    % Storing d_k avoids an extra subtraction to obtain s_{k-1} in the next
    % iteration
    parameters.dk = dk;
else
    % Following iterations
    % Load variables
    yk = parameters.yk;
    dk = parameters.dk;
    gk = parameters.gk;
    
    % Compute alpha_k
    if isfield(parameters,'type') && parameters.type == 2
        alpha = (yk' * dk) / (yk' * yk);
    else
        alpha = (dk' * dk) / (dk' * yk);
    end
    
    
    % New updated step
    dk = -alpha * gk;
    
    % Barzilai-Borwein gradient descent iteration
    x = x + dk;
    
    % Store next gradient
    parameters.gk = grad(x);
    parameters.yk = parameters.gk - gk;
    parameters.dk = dk;    
end



end

