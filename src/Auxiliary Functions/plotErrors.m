function h = plotErrors(errors, names, yaxis)
%PLOTERRORS Produces plots of errors with y-axis in logarithm scale.
%   Accepts a cell array of all errors of each of the algorithm specified
%   in the vector names and plots all of them.

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

h = figure;
for currentAlgorithm = 1:length(errors)
    semilogy(errors{currentAlgorithm});
    hold on;
end
legend(names);
xlabel('iterations');
ylabel(yaxis,'Interpreter','latex');

