function w_max = optDOR(A,functionT)
%OPTDOR Returns the optimal w that guarantees the fastest convergence for
%the Delayed Over-relaxation method. 
%   Solves the optimization problem:
%       min    max{ lambda( T_DOR ) }
%   Given the format of the problem, the solution is found through
%   bissection.

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

% % ==== Brute Force Option ====
% 
% G = functionT(A);
% 
% N = 200;
% w_values = linspace(1,1.3,N);
% spectral_radius = zeros(1,N);
% 
% rho_min = 1;
% 
% for i = 1:N
%     w = w_values(i);
%     T = [w*G, (1-w)*eye(length(A));
%          eye(length(A)), zeros(length(A))];
%     
%     lambdaT = eig(T);
%     
%     spectral_radius(i) = max(abs(lambdaT));    
%     
%     if spectral_radius(i) < rho_min
%         rho_min = spectral_radius(i);
%         w_max = w;
%     end
% end
% 
% plot(w_values,spectral_radius);


% ==== Bissection ====
w_max = optQuasiConvex(A, @Tdor, functionT);
end

