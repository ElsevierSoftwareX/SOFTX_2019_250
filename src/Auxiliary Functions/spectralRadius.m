function r = spectralRadius(A)
%SPECTRALRADIUS Computes the spectral radius of a matrix A.
%   Calculates r = \rho(A) = max |\lambda(A)|.

r = max(abs(eig(A))); 
end

