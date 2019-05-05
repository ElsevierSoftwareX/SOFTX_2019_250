function [algorithms, parameters, names] = getParameters(A, b, r, methods)
%GETPARAMETERS Returns the data structures with (optimal when possible)
%parameters for the various methods specified in 'methods'. If called
%with no argument it prints all possible algorithm names.
%   Returns data structures for optimization algorithms solving
%       f(x) = 0.5 x'Qx - p' x + r
%   where Q = A'A and p = A'b.
%   Additionally, for linear equation solvers applied to
%       grad f(x) = 0 <=> Qx = p
%   for algorithms requiring symmetric matrices or
%       Ax = b.
%   The parameters for the optimization algorithms follow the idea in [14].

% Add descriptions for new implementations in the toolbox
%                 Names             parameter names
optAlgorithmNames = ["Gradient Descent","alpha"
    "Heavy-Ball","alpha, beta"
    "Nesterov","alpha, beta"
    "Delayed Nesterov","alpha, beta"
    "Momentum","alpha, beta"
    "FISTA","alpha"
    "DFISTA","alpha, f(x)"
    "Second Nesterov","alpha"
    "Barzilai-Borwein","no parameters"
    "Random Descent","no parameters"
    "Cauchy-Barzilai-Borwein","Q"
    "General Barzilai-Borwein","f(x), M, rho, delta, sigma, amin, amax"];

linEqAlgorithmNames = ["Jacobi","no parameters"
    "WeightedJacobi","w"
    "Gauss-Seidel","no parameters"
    "Successive Over-relaxation","w"
    "Richardson","w"
    "Conjugate Gradient","no parameters"
    "Biconjugate Gradient","no parameters"
    "Newton-Raphson","no parameters"
    "Sparse Broyden","no parameters"
    "Broyden","no parameters"
    "Bad Broyden", "no parameters"
    "Delayed Over-Relaxation", "method, w"
    "Minimal Residual - DOR","no parameters"
    "Accelerated Over-Relaxation","w, r"
    "PAOSOR","no parameters"
    "Alternating Anderson-Jacobi","m, p, beta, w"
    "Chebyshev","no parameters"
    "Quasi-Chebyshev","alpha, jacobi"
    "HSS","alpha"
    "Kaczmarz","no parameters"
    "Coordinate Descent","no parameters"
    "CGNE","no parameters"
    "IBiCG","no parameters"];

% Optimization Algorithms
switch nargin
    case 0
        % If called with no input arguments display help
        fprintf('Select \''all\'' for all implemented algorithms.\n The available optimization algorithms and their parameters are:\n');
        fprintf('\t %30s \t with %s\n', optAlgorithmNames');
        fprintf('The available optimization algorithms and their parameters are:\n');
        fprintf('\t %30s \t with %s\n', linEqAlgorithmNames');
    otherwise
        % Are we selecting all methods
        allMethods = any(strcmp(methods,'all'));
        
        if allMethods
            % Make sure there is only one element in methods
            methods = 'a';
            
            % Already declare the names variable
            names.optimization = optAlgorithmNames(:,1)';
            names.linearEquations = linEqAlgorithmNames(:,1)';
            
            optAlgorithms = cell(1,length(optAlgorithmNames));
            optParameters = cell(1,length(optAlgorithmNames));
            
            linEqAlgorithms = cell(1,length(linEqAlgorithmNames));
            linEqParameters = cell(1,length(linEqAlgorithmNames));
        else
            typeOptimization = ismember(methods, optAlgorithmNames);
            typeLinearSolver = ismember(methods, linEqAlgorithmNames);
            
            names.optimization = strings(1,sum(typeOptimization));
            names.linearEquations = strings(1,sum(typeLinearSolver));
            
            optAlgorithms = cell(1,length(intersect(methods,optAlgorithmNames)));
            optParameters = cell(1,length(intersect(methods,optAlgorithmNames)));
            
            linEqAlgorithms = cell(1,length(intersect(methods,linEqAlgorithmNames)));
            linEqParameters = cell(1,length(intersect(methods,linEqAlgorithmNames)));
        end
        
        % indices to write on data structures
        optIndex = 1;
        linEqIndex = 1;
        
        % Important quantities regarding the Optimization problem at hand
        Q = A'*A;
        lambdaQ = eig(Q);
        % Strong convexity parameter
        m = min( lambdaQ );
        % Lipschitz gradient parameter
        L = max( lambdaQ );
        % Condition ration of the function
        kappa = L/m;
        
        % Important quantities regarding the Linear Equation problem
        lambda_min = min(eig(A));
        lambda_max = max(eig(A));
        
        % Let us build the data structures with optimal parameters when
        % known for each of the elements in methods
        for i = 1:length(methods)
            %==== Test if name is an Optimization Algorithm ====
            if ~allMethods && typeOptimization(i)
                names.optimization(optIndex) = methods(i);
            %==== Test if name is an Optimization Algorithm ====
            elseif ~allMethods && typeLinearSolver(i)
                names.linearEquations(linEqIndex) = methods(i);
            end
            
            %==== Optimization Algorithms ====
            % Gradient Descent
            if strcmp(methods(i),"Gradient Descent") || allMethods
                optAlgorithms{optIndex} = @gradientDescent;
                optParameters{optIndex} = struct('alpha',2/(L + m));
                optIndex = optIndex + 1;
            end
            % Heavy-Ball
            if strcmp(methods(i),"Heavy-Ball") || allMethods
                optAlgorithms{optIndex} = @heavyBall;
                optParameters{optIndex} = struct('alpha',4/(sqrt(L) + sqrt(m))^2,'beta',( (sqrt(kappa) - 1)/(sqrt(kappa) + 1) )^2);
                optIndex = optIndex + 1;
            end
            % Nesterov
            if strcmp(methods(i),"Nesterov") || allMethods
                optAlgorithms{optIndex} = @nesterov;
                optParameters{optIndex} = struct('alpha',4/(3*L + m),'beta',(sqrt(3*kappa + 1) - 2)/(sqrt(3*kappa + 1) + 2) );
                optIndex = optIndex + 1;
            end
            % Delayed Nesterov (Not Optimal)
            if strcmp(methods(i),"Delayed Nesterov") || allMethods
                optAlgorithms{optIndex} = @delayedNesterov;
                optParameters{optIndex} = struct('alpha',4/(3*L + m),'beta',(sqrt(3*kappa + 1) - 2)/(sqrt(3*kappa + 1) + 2) );
                optIndex = optIndex + 1;
            end
            % Momentum (Not Optimal)
            if strcmp(methods(i),"Momentum") || allMethods
                optAlgorithms{optIndex} = @momentum;
                optParameters{optIndex} = struct('alpha',4/(3*L + m),'beta',(sqrt(3*kappa + 1) - 2)/(sqrt(3*kappa + 1) + 2) );
                optIndex = optIndex + 1;
            end
            % Fast Iterative Shrinkage-Thresholding Algorithm (Not Optimal)
            if strcmp(methods(i),"FISTA") || allMethods
                optAlgorithms{optIndex} = @fista;
                optParameters{optIndex} = struct('alpha',4/(3*L + m));
                optIndex = optIndex + 1;
            end
            % Descent Fast Iterative Shrinkage-Thresholding Algorithm (Not Optimal)
            if strcmp(methods(i),"DFISTA") || allMethods
                optAlgorithms{optIndex} = @dfista;
                optParameters{optIndex} = struct('alpha',4/(3*L + m),'f',@(x) 0.5*x'*Q*x - b'*A*x + r);
                optIndex = optIndex + 1;
            end
            % Second Nesterov (Not Optimal)
            if strcmp(methods(i),"Second Nesterov") || allMethods
                optAlgorithms{optIndex} = @secondNesterov;
                optParameters{optIndex} = struct('alpha',4/(3*L + m));
                optIndex = optIndex + 1;
            end
            % Barzilai-Borwein (Could Optimize first alpha)
            if strcmp(methods(i),"Barzilai-Borwein") || allMethods
                optAlgorithms{optIndex} = @barzilaiBorwein;
                optParameters{optIndex} = [];
                optIndex = optIndex + 1;
            end
            % Random Descent (Not Optimal)
            if strcmp(methods(i),"Random Descent") || allMethods
                optAlgorithms{optIndex} = @randomDescent;
                optParameters{optIndex} = [];
                optIndex = optIndex + 1;
            end
            % Cauchy-Barzilai-Borwein
            if strcmp(methods(i),"Cauchy-Barzilai-Borwein") || allMethods
                optAlgorithms{optIndex} = @cbb;
                optParameters{optIndex} = struct('Q',Q);
                optIndex = optIndex + 1;
            end
            % General Barzilai-Borwein (Not Optimal)
            if strcmp(methods(i),"General Barzilai-Borwein") || allMethods
                optAlgorithms{optIndex} = @gbb;
                optParameters{optIndex} = struct('f',@(x) 0.5*x'*Q*x - b'*A*x + r,'M',10,'rho',0.5,'delta',0.3,'sigma',0.8,'amin',0.1,'amax',10);
                optIndex = optIndex + 1;
            end
            %==== Linear equations Solvers ====
            % Jacobi
            if strcmp(methods(i),"Jacobi") || allMethods
                linEqAlgorithms{linEqIndex} = @jacobi;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % WeightedJacobi
            if strcmp(methods(i),"WeightedJacobi") || allMethods
                linEqAlgorithms{linEqIndex} = @weightedJacobi;
                linEqParameters{linEqIndex} = struct('w',optWJ(lambda_min,lambda_max));
                linEqIndex = linEqIndex + 1;
            end
            % Gauss-Seidel
            if strcmp(methods(i),"Gauss-Seidel") || allMethods
                linEqAlgorithms{linEqIndex} = @gaussSeidel;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Successive Over-relaxation
            if strcmp(methods(i),"Successive Over-relaxation") || allMethods
                linEqAlgorithms{linEqIndex} = @sor;
                linEqParameters{linEqIndex} = struct('w',optSOR(A));
                linEqIndex = linEqIndex + 1;
            end
            % Richardson
            if strcmp(methods(i),"Richardson") || allMethods
                linEqAlgorithms{linEqIndex} = @richardson;
                linEqParameters{linEqIndex} = struct('w',optRichardson(lambda_min,lambda_max));
                linEqIndex = linEqIndex + 1;
            end
            % Conjugate Gradient
            if strcmp(methods(i),"Conjugate Gradient") || allMethods
                linEqAlgorithms{linEqIndex} = @conjugateGradient;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Biconjugate Gradient
            if strcmp(methods(i),"Biconjugate Gradient") || allMethods
                linEqAlgorithms{linEqIndex} = @biconjugateGradient;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Newton-Raphson
            if strcmp(methods(i),"Newton-Raphson") || allMethods
                linEqAlgorithms{linEqIndex} = @newton;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Sparse Broyden
            if strcmp(methods(i),"Sparse Broyden") || allMethods
                linEqAlgorithms{linEqIndex} = @sparseBroyden;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Broyden
            if strcmp(methods(i),"Broyden") || allMethods
                linEqAlgorithms{linEqIndex} = @broyden;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Bad Broyden
            if strcmp(methods(i),"Bad Broyden") || allMethods
                linEqAlgorithms{linEqIndex} = @badBroyden;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Delayed Over-Relaxation (Not Working)
            if strcmp(methods(i),"Delayed Over-Relaxation") || allMethods
                linEqAlgorithms{linEqIndex} = @dor;
                linEqParameters{linEqIndex} = struct('method',@jacobi,'parameters',[],'w',optDOR(A));
                linEqIndex = linEqIndex + 1;
            end
            % Minimal Residual - DOR
            if strcmp(methods(i),"Minimal Residual - DOR") || allMethods
                linEqAlgorithms{linEqIndex} = @mrDOR;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Accelerated Over-Relaxation
            if strcmp(methods(i),"Accelerated Over-Relaxation") || allMethods
                linEqAlgorithms{linEqIndex} = @aor;
                [w, r] = optAOR(A);
                linEqParameters{linEqIndex} = struct('w',w,'r',r);
                linEqIndex = linEqIndex + 1;
            end
            % PAOSOR
            if strcmp(methods(i),"PAOSOR") || allMethods
                linEqAlgorithms{linEqIndex} = @paosor;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Alternating Anderson-Jacobi (Not Optimal)
            if strcmp(methods(i),"Alternating Anderson-Jacobi") || allMethods
                linEqAlgorithms{linEqIndex} = @aaj;
                linEqParameters{linEqIndex} = struct('m',10,'p',6,'beta',0.3,'w',optWJ(lambda_min,lambda_max));
                linEqIndex = linEqIndex + 1;
            end
            % Chebyshev
            if strcmp(methods(i),"Chebyshev") || allMethods
                linEqAlgorithms{linEqIndex} = @chebyshev;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Quasi-Chebyshev
            if strcmp(methods(i),"Quasi-Chebyshev") || allMethods
                linEqAlgorithms{linEqIndex} = @quasiChebyshev;
                linEqParameters{linEqIndex} = struct('alpha',1,'jacobi',1);
                linEqIndex = linEqIndex + 1;
            end
            % HSS
            if strcmp(methods(i),"HSS") || allMethods
                linEqAlgorithms{linEqIndex} = @hss;
                linEqParameters{linEqIndex} = struct('alpha',optHSS(A,1));
                linEqIndex = linEqIndex + 1;
            end
            % Kaczmarz
            if strcmp(methods(i),"Kaczmarz") || allMethods
                linEqAlgorithms{linEqIndex} = @kaczmarz;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % Coordinate Descent
            if strcmp(methods(i),"Coordinate Descent") || allMethods
                linEqAlgorithms{linEqIndex} = @coordinateDescent;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % CGNE
            if strcmp(methods(i),"CGNE") || allMethods
                linEqAlgorithms{linEqIndex} = @cgne;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
            % IBiCG
            if strcmp(methods(i),"IBiCG") || allMethods
                linEqAlgorithms{linEqIndex} = @improvedBiCG;
                linEqParameters{linEqIndex} = [];
                linEqIndex = linEqIndex + 1;
            end
        end
        
        % Store output
        algorithms.optimization = optAlgorithms;
        algorithms.linearEquations = linEqAlgorithms;
        
        parameters.optimization = optParameters;
        parameters.linearEquations = linEqParameters;
        
end

