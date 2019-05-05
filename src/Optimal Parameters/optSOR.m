function w_max = optSOR(A)
%OPTSOR Returns the optimal w that guarantees the fastest convergence for
%the Successive Over-relaxation method. 
%   Solves the optimization problem:
%       min    max{ lambda( T_SOR ) }
%   Given the format of the problem, the solution is found through brute
%   force.

% ==== Could NOT solve using optimization (to research) ====
% ==== Brute Force Option ====
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

N = 100;
w_values = linspace(0,2,N);
spectral_radius = zeros(1,N);

rho_min = 1;
w_max = 1;

for i = 1:N
    w = w_values(i);
    T = -(D-D*(w*L)^1+D*(w*L)^2-D*(w*L)^3)*(w*U + (w-1)*D);
    
    lambdaT = eig(T);
    
    spectral_radius(i) = max(abs(lambdaT));    
    
    if spectral_radius(i) < rho_min
        rho_min = spectral_radius(i);
        w_max = w;
    end
end

% plot(1:N,spectral_radius);
% ==== Closed-form solution ====
% 
end

