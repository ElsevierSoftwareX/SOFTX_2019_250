function C = Csor(A,b,w)
%CSOR Computes the vector in the matrix form of the SOR algorithm.
%   Computes :
%   C = (D+(w*L))^-1 * wb 

D = diag(diag(A));
L = tril(A,-1);

C = (D+w*L)\w*b;

end

