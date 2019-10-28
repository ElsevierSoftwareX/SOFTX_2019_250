function [ stateVectors , errors ] = linSolver( algorithms , parameters , A , b , errorFunction, max_iterations , initialState, projectionFunction, tol)
%LINSOLVER Applies the iterations of all 'algorithms' using their respective 'parameters'.
%   Applies all algorithms using their parameters and produces the state
%   evolutions and their respective errors given the provided error function.

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

% Check whether there are parameters for all algorithms
if length(algorithms) ~= length(parameters)
    ME = MException('linSolver:mismatchParameters',...
                        'The size of the list of Algorithms and Parameters is different');
    throw(ME);
end

% If no tolerance was given set it to 1E-14
if ~exist('tol','var')
    tol = 1E-14;
end

% Initialize the Outputs
stateVectors = cell(1,length(algorithms));
errors = cell(1,length(algorithms));


% Cycle over all the specified algorithms
for currentAlgorithm = 1:length(algorithms)
    
    % State and error vectors
    x = zeros(size(initialState,1),max_iterations);
    error = zeros(1,max_iterations);
    
    % Initialize state and error measure
    x(:,1) = initialState;
    error(1) = errorFunction(initialState);
    
    % Simulation of all allowed iterations
    for time = 2:max_iterations
        
        % Execute one step of the current algorithm
        if time == 2 
            % It is the first iteration and there is no x(:,time-2)
            [x(:,time) , parameters{currentAlgorithm}] = algorithms{currentAlgorithm}(x(:,time-1) , x(:,time-1) , A , b , parameters{currentAlgorithm} , max_iterations);
        else
            [x(:,time) , parameters{currentAlgorithm}] = algorithms{currentAlgorithm}(x(:,time-1) , x(:,time-2) , A , b , parameters{currentAlgorithm} , max_iterations);
        end
        
        % If a specified projection function is supplied the state is
        % projected
        if exist('projectionFunction','var')
            x(:,time) = projectionFunction( x(:,time) );
        end
        
        % Compute the new error of the current state
        error(time) = errorFunction(x(:,time));
        
        % If given tolerance is met then go to the next algorithm
        if error(time) <= tol
            error(time+1:end) = repmat(error(time),[1 max_iterations-time]);
            x(:,time+1:end) = repmat(x(:,time),[1 max_iterations-time]);
            break;
        end
        
    end
    
    % Add the state and error vectors to the output
    stateVectors{currentAlgorithm} = x;
    errors{currentAlgorithm} = error;
    
end
end

