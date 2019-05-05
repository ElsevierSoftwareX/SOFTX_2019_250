clearvars;
close all;


% Matrix used in this example
% A = [0 0 0 1/3;
%     1 0 1/2 1/3;
%     0 1/2 0 1/3;
%     0 1/2 1/2 0];

A = adjacencyBA(20,5,5);

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
optAlgorithmNames = ["Nesterov","Second Nesterov","FISTA","Barzilai-Borwein","Random Descent"];
% From Linear Equation Solvers
% linEqAlgorithmNames = ["Jacobi","WeightedJacobi","Gauss-Seidel","Successive Over-relaxation","Richardson","Conjugate Gradient","Biconjugate Gradient","Newton-Raphson","Sparse Broyden","Broyden","Bad Broyden", "Delayed Over-Relaxation", "Minimal Residual - DOR","Accelerated Over-Relaxation","PAOSOR","Alternating Anderson-Jacobi","Chebyshev","Quasi-Chebyshev","HSS","Kaczmarz","Coordinate Descent","CGNE","IBiCG"];
linEqAlgorithmNames = ["Jacobi","Successive Over-relaxation","WeightedJacobi","Gauss-Seidel"];
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