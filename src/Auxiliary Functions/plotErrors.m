function h = plotErrors(errors, names, yaxis)
%PLOTERRORS Produces plots of errors with y-axis in logarithm scale.
%   Accepts a cell array of all errors of each of the algorithm specified
%   in the vector names and plots all of them.

h = figure;
for currentAlgorithm = 1:length(errors)
    semilogy(errors{currentAlgorithm});
    hold on;
end
legend(names);
xlabel('iterations');
ylabel(yaxis,'Interpreter','latex');

