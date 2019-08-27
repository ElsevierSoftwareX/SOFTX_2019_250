clearvars;
close all;

A = adjacencyBA(5,5,5);

% Typical parameter values used in defining the PageRank problem
m = 0.15;

% Parameters for the simulation
networkSizes = 100:100:200;
stoppingTimes = zeros(3,length(networkSizes));
tolerance = 1E-8;
optimalW = zeros(1,length(networkSizes));
chosenW = zeros(1,length(networkSizes));

for epoch = 1:length(networkSizes)
    n = networkSizes(epoch);
    A = growAdjacencyBA(n,double(logical(A)),5);
    
    M = (1-m)*A + (m/n)*ones(n);
    
    % Centralized Solution - solved by finding the eigenvector
    [V, D] = eig(M);
    [lambda_max, I] = max(diag(D));
    
    % Solution
    x_star = abs(V(:,I)/norm(V(:,I),1));
    
    % Number of maximum iterations for the simulation
    max_iterations = 100;
    
    % Define linear gradient function
    gradA = eye(n)-(1-m)*A;
    gradb = (m/n)*ones(n,1);
    
    % ------- Algorithms names to be tested --------
    linEqAlgorithmNames = ["Jacobi","WeightedJacobi","Successive Over-relaxation"];
    % Get parameters and function handles
    [algorithms, parameters, names] = getParameters(gradA,gradb,0.5*(gradb'*gradb),linEqAlgorithmNames);
    
    % Update parameter w for SOR
    optimalW(epoch) = parameters.linearEquations{3}.w;
    if epoch >2
        parameters.linearEquations{3}.w = chosenW(1,epoch-1) + 1E-2;
    end
    chosenW(epoch) = parameters.linearEquations{3}.w;
    % Define an error functions
    errorFunction = @(x) norm(x-x_star,2);
    
    % Use the tool for linear equation solvers to run all algorithms
    [ linstateVectors , linerrors ] = linSolver(algorithms.linearEquations , parameters.linearEquations , gradA , gradb , errorFunction , max_iterations , ones(n,1)/n, @(x) x, tolerance);
    
    for i = 1:length(linerrors)
        stopped = find(linerrors{i}<=1E-8,1);
        if isempty (stopped)
            stopped = max_iterations;
        end
        stoppingTimes(i,epoch) = stopped;
    end
end

figure;
h = plot(networkSizes,stoppingTimes(1,:),'-o','LineWidth',2);
hold on;
plot(networkSizes,stoppingTimes(2,:),'-s','LineWidth',2);
hold on;
plot(networkSizes,stoppingTimes(3,:),'-*','LineWidth',2);
legend(names.linearEquations);
xlabel('number of pages');
ylabel('Stopping iteration');
% saveas(h, 'growingBAnetLinearW','pdf');

figure;
plot(optimalW);
hold on;
plot(chosenW);

% Plot the various error functions for the linear equation solver algorithms
figure;
for currentAlgorithm = 1:length(algorithms.linearEquations)
    semilogy(linerrors{currentAlgorithm});
    hold on;
end
legend(names.linearEquations);
xlabel('iterations');
ylabel('$\|x(k) - x_\infty \|_2$','Interpreter','latex');
