function [ stateVectors , errors ] = quadSolver( algorithmNames , parameters , A , b , errorFunction, max_iterations , initialState ,  experimentName, projectionFunction, tol, errorDescription)
%QUADSOLVER Applies the iterations of all optimization and linear equation
%   solvers to the quadratic problem. 
%
%   It solves:
%   0.5 x'A'A x - (A'b)' x + 0.5 b'b 
%   which is equivalent to;
%   0.5 ||Ax - b||^2 
%   It applies the optimization algorithms and linear equation solvers both
%   for Ax = b and the normal equation A'Ax = A'b.

% Get the algorithms functions and optimal parameters for the standard
% version
[algorithms, optimalParameters, names] = getParameters(A,b,0.5*(b'*b),algorithmNames);

% If no tolerance was given set it to 1E-14
if ~exist('tol','var')
    tol = 1E-14;
end

% If no projection is specified we can define the identity function
if ~exist('projectionFunction','var')
    projectionFunction = @(x) x;
end

% ====== Optimization Algorithms ======
% Gradient function of the quadratic function
% Define gradient function
grad = @(x) A'* (A*x - b);
if isfield(parameters,'optimization') && length(parameters.optimization) == length(optimalParameters.optimization)
    % The parameters from the user are going to be used
    [ stateVectors.optimization , errors.optimization ] = optSolver(algorithms.optimization , parameters.optimization , grad , errorFunction , max_iterations , initialState , projectionFunction , tol);
else
    fprintf("Optimization parameter not correctly specified. Selected optimal ones.\n");
    [ stateVectors.optimization , errors.optimization ] = optSolver(algorithms.optimization , optimalParameters.optimization , grad , errorFunction , max_iterations , initialState , projectionFunction , tol);
end

% ====== Linear Equations Solvers ======
if isfield(parameters,'linearEquations') && length(parameters.linearEquations) == length(optimalParameters.linearEquations)
    % The parameters from the user are going to be used
    [ stateVectors.linearEquations , errors.linearEquations ] = optSolver(algorithms.linearEquations , parameters.linearEquations , A , b , errorFunction , max_iterations , initialState , projectionFunction , tol);
else
    fprintf("Equation Solvers parameter not correctly specified. Selected optimal ones.\n");
    [ stateVectors.linearEquations , errors.linearEquations ] = linSolver(algorithms.linearEquations , optimalParameters.linearEquations , A , b , errorFunction , max_iterations , initialState , projectionFunction , tol);
end

% ====== Normal Equations Solvers ======
if isfield(parameters,'linearEquations') && length(parameters.linearEquations) == length(optimalParameters.linearEquations)
    % The parameters from the user are going to be used
    [ stateVectors.normalEquations , errors.normalEquations ] = optSolver(algorithms.linearEquations , parameters.linearEquations , A'*A , A'*b , errorFunction , max_iterations , initialState , projectionFunction , tol);
else
    fprintf("Normal Equation Solvers parameter not correctly specified. Selected optimal ones.\n");
    [ stateVectors.normalEquations , errors.normalEquations ] = linSolver(algorithms.linearEquations , optimalParameters.linearEquations , A'*A , A'*b , errorFunction , max_iterations , initialState , projectionFunction , tol);
end

% ====== Plot and Save the Results ======
if exist('experimentName','var')
    experimentName = strcat(userpath, '/OPTool/Stored Outputs/', experimentName);
    % If no function description was provided convert the errorFunction to
    % string
    if ~exist('errorDescription','var')
        errorDescription = strcat("$",strrep(char(errorFunction),'@(x)',''),"$");
    end
    % Plot optimization, linear equation and normal equation
    optFigure = plotErrors(errors.optimization , names.optimization , errorDescription);
    linEqFigure = plotErrors(errors.linearEquations,names.linearEquations,errorDescription);
    normalFigure = plotErrors(errors.normalEquations,names.linearEquations,errorDescription);
    % Store the variables
    save(experimentName,'stateVectors','errors');
    % Save figures
    saveas(optFigure, strcat(experimentName, '-optimization'),'pdf');
    saveas(linEqFigure,strcat(experimentName, '-linEqSolvers'),'pdf');
    saveas(normalFigure,strcat(experimentName, '-normalEquations'),'pdf');
end
