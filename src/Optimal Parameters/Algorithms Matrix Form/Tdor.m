function T = Tdor(A,G,w)
%TDOR Computes the transition matrix Tsor for the DOR algorithm.
%   Computes :
%   T = [ wG (1-w)I; I 0]

n = length(A);

T = [ w*G(A) (1-w)*eye(n); eye(n) zeros(n)];

end

