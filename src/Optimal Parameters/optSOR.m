function w_max = optSOR(A)
%OPTSOR Returns the optimal w that guarantees the fastest convergence for
%the Successive Over-relaxation method. 
%   Solves the optimization problem:
%       min    max{ lambda( T_SOR ) }
%   Given the format of the problem, the solution is found through
%   bissection.

% % ==== Brute Force Option ====
% D = diag(diag(A));
% L = tril(A,-1);
% U = triu(A,1);
% 
% N = 400;
% w_values = linspace(0,2,N);
% spectral_radius = zeros(1,N);
% 
% rho_min = 1;
% w_max = 1;
% 
% for i = 1:N
%     w = w_values(i);
%     
%     % The transition matrix T is given by:
%     % T = -(D+(w*L))^-1*(w*U + (w-1)*D) 
%     % where we explore the fact that:
%     % iLw = (D+(w*L))^-1 = \sum_{j=0}^{n-1} (-1)^j D*(w*L)^j
% %     iLw = D;
% %     for j = 1:length(A)-1
% %         iLw = iLw + (-1)^j* D*(w*L)^j;
% %     end
% %     T = -iLw*(w*U + (w-1)*D);
% 
%     T = -inv(D+(w*L))*(w*U + (w-1)*D);
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
% h = figure;
% plot(w_values,spectral_radius);
% xlabel('$\omega$','Interpreter','latex');
% ylabel('$\rho(T_{sor}(\omega))$','Interpreter','latex');
% saveas(h, 'quasiConvexSOR','pdf');

% ==== Bissection ====
w_max = optQuasiConvex(A,@Tsor,[]);

end

