clearvars;
close all;

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


% Matrix used in this example
% A = [0 0 0 1/3;
%     1 0 1/2 1/3;
%     0 1/2 0 1/3;
%     0 1/2 1/2 0];

A = adjacencyBA(200,5,5);

% Typical parameter values used in defining the PageRank problem
m = 0.15;
n = size(A,1);

M = (1-m)*A + (m/n)*ones(n);

% Centralized Solution - solved by finding the eigenvector 
[V, D] = eig(M);
[lambda_max, I] = max(diag(D));

% Solution
x_star = abs(V(:,I)/norm(V(:,I),1));

% Number of maximum iterations for the simulation
max_iterations = 50;

% Define gradient function
grad = @(x) ((1-m)*A-eye(n))'*(((1-m)*A-eye(n))*x+(m/n)*ones(n,1));

% Define linear gradient function
gradA = eye(n)-(1-m)*A;
gradb = (m/n)*ones(n,1);

% ------- Algorithms names to be tested --------
% From Optimization
% optAlgorithmNames = ["Gradient Descent", "Heavy-Ball", "Nesterov","Delayed Nesterov","Momentum","FISTA","DFISTA","Second Nesterov","Barzilai-Borwein","Random Descent","Cauchy-Barzilai-Borwein","General Barzilai-Borwein"];
optAlgorithmNames = ["Gradient Descent", "Nesterov","Barzilai-Borwein","Cauchy-Barzilai-Borwein"];
% From Linear Equation Solvers
% linEqAlgorithmNames = ["Jacobi","WeightedJacobi","Gauss-Seidel","Successive Over-relaxation","Richardson","Conjugate Gradient","Biconjugate Gradient","Newton-Raphson","Sparse Broyden","Broyden","Bad Broyden", "Delayed Over-Relaxation", "Minimal Residual - DOR","Accelerated Over-Relaxation","PAOSOR","Alternating Anderson-Jacobi","Chebyshev","Quasi-Chebyshev","HSS","Kaczmarz","Coordinate Descent","CGNE","IBiCG"];
linEqAlgorithmNames = ["Jacobi","WeightedJacobi","Successive Over-relaxation"];
% Get parameters and function handles
[algorithms, parameters, names] = getParameters(gradA,gradb,0.5*(gradb'*gradb),[optAlgorithmNames linEqAlgorithmNames]);

% Define an error functions
errorFunction = @(x) norm(x-x_star,2);

% Use the tool for optimization solvers to run all algorithms
[ stateVectors , errors ] = optSolver(algorithms.optimization , parameters.optimization , grad , errorFunction , max_iterations , ones(n,1)/n);
% Use the tool for linear equation solvers to run all algorithms
[ linstateVectors , linerrors ] = linSolver(algorithms.linearEquations , parameters.linearEquations , gradA , gradb , errorFunction , max_iterations , ones(n,1)/n);

% Plot the various error functions for the optimization algorithms
figure;
for currentAlgorithm = 1:length(algorithms.optimization)
    semilogy(errors{currentAlgorithm});
    hold on;
end
legend(names.optimization);
xlabel('iterations');
ylabel('$\|x(k) - x_\infty \|_2$','Interpreter','latex');

% Plot the various error functions for the linear equation solver algorithms
figure;
for currentAlgorithm = 1:length(algorithms.linearEquations)
    semilogy(linerrors{currentAlgorithm});
    hold on;
end
legend(names.linearEquations);
xlabel('iterations');
ylabel('$\|x(k) - x_\infty \|_2$','Interpreter','latex');

% ==== Plot Optimization vs Linear Equation ==== %
h = figure;
for currentAlgorithm = 1:length(algorithms.optimization)
    semilogy(errors{currentAlgorithm});
    hold on;
end
semilogy(linerrors{1});
legend([names.optimization names.linearEquations(1)], 'Location', 'southwest');
xlabel('iterations');
ylabel('$\|x(k) - x_\infty \|_2$','Interpreter','latex');

% saveas(h, 'comparisonJacobiOptimization','pdf');