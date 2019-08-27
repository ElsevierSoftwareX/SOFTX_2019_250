function A = growAdjacencyBA(N, S, m)
%ADJACENCYBA Generates a random scale-free Markov Chain matrix using
%Barabási–Albert model.
%   The random network has nodes being added with higher probability to
%   high-degree nodes already in the network.

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

