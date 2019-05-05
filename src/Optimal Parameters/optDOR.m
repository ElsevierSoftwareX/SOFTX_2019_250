function w_max = optDOR(A)
%OPTDOR Returns the optimal w that guarantees the fastest convergence for
%the Delayed Over-relaxation method. 
%   Solves the optimization problem:
%       min    

% ==== Could NOT solve using optimization (to research) ====
% ==== Brute Force Option ====
D = diag(diag(A));
R = A - D;

G = -inv(D)*R;

N = 100;
w_values = linspace(0,2,N);
spectral_radius = zeros(1,N);

rho_min = 1;

for i = 1:N
    w = w_values(i);
    T = [w*G, (1-w)*eye(length(A));
         eye(length(A)), zeros(length(A))];
    
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

