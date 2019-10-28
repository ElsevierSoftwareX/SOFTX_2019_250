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


% Number of nodes in the network
n = 5;

% Matrix D and vectors defining the g(x) = 0.5||Dx-v 1_n+e_n||^2
D = -eye(n) + diag(ones(1,n-1),1) + e(n,n)*e(1,n)';
S = -eye(n) + diag(ones(1,n-1),1);
v = 1/n;

en = e(n,n);

% Number of maximum iterations for the simulation
% max_iterations = round(n/10)*50;
max_iterations = 20;

% Define gradient function
grad = @(x) (D'*D*x+D'*en);

% Define linear gradient function
gradA = D;
gradb = ones(n,1)/n-en;

% ------- Algorithms names to be tested --------
% From Optimization
optAlgorithmNames = ["Gradient Descent", "Nesterov","NesterovLTV","Heavy-Ball"];
% From Linear Equation Solvers
linEqAlgorithmNames = "Gauss-Seidel";
% Get parameters and function handles
[algorithms, parameters, names] = getParameters(gradA,gradb,0.5*(gradb'*gradb),[optAlgorithmNames linEqAlgorithmNames]);

% Define an error functions
errorFunction = @(x) norm(S*(x-min(x)*ones(n,1))+en-ones(n,1)/n);

% Use the tool for optimization solvers to run all algorithms
[ stateVectors , errors ] = optSolver(algorithms.optimization , parameters.optimization , grad , errorFunction , max_iterations , ones(n,1)/n);
% Use the tool for linear equation solvers to run all algorithms
[ linstateVectors , linerrors ] = linSolver(algorithms.linearEquations , parameters.linearEquations , D'*D , -D'*en , errorFunction , max_iterations , ones(n,1)/n);

% Plot the various error functions for the optimization algorithms
h = figure;
for currentAlgorithm = 1:length(algorithms.optimization)
    semilogy(errors{currentAlgorithm});
    hold on;
end
semilogy(linerrors{1})
hlegend = legend([names.optimization "Gauss-Seidel"]);
xlabel('iterations');
ylabel('$\|\phi(k) - \phi^\star\|$','Interpreter','latex');
set(hlegend,'Location','best');
saveas(h,strcat("GSvsOptimization", num2str(n)),'pdf');
