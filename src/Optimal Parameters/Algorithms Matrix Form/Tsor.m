function T = Tsor(A,~,w)
%TSOR Computes the transition matrix Tsor for the SOR algorithm.
%   Computes :
%   T = -(D+(w*L))^-1*(w*U + (w-1)*D) 
% where we explore the fact that:
%   iLw = (D+(w*L))^-1 = \sum_{j=0}^{n-1} (-1)^j D*(w*L)^j

D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);

% iLw = D;
% for j = 1:length(A)-1
%     iLw = iLw + (-1)^j* D*(w*L)^j;
% end
% T = -iLw*(w*U + (w-1)*D);

T = -inv(D+(w*L))*(w*U + (w-1)*D);

end

