function C = Cj(A,b)
%CJ Computes the vector in the matrix form of the Jacobi algorithm.
%   Computes :
%   C = D^-1*b 

D = diag(diag(A));

C = D\b;

end

