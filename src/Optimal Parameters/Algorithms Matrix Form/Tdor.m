function T = Tdor(A,G,w)
%TDOR Computes the transition matrix Tsor for the DOR algorithm.
%   Computes :
%   T = [ wG (1-w)I; I 0]

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

n = length(A);

T = [ w*G(A) (1-w)*eye(n); eye(n) zeros(n)];

end

