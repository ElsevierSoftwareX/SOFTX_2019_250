function A = growAdjacencyBA(N, S, m)
%ADJACENCYBA Generates a random scale-free Markov Chain matrix using
%Barabási–Albert model.
%   The random network has nodes being added with higher probability to
%   high-degree nodes already in the network.

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

% Generate a random BA topology
network = growBAgraph(N, S, m);

% Count the nodes degree
degree = sum(network);

% Check for dangling nodes
index = (degree == 0);
% Account for dangling nodes
network(:,index) = network(index,:)';

% Compute the Markov Chain matrix
A = network./sum(network);
end

