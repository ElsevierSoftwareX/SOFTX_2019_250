function CNet(Net)

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

format compact
format long e
theta = linspace(0,2*pi,length(Net)+1);
xy = zeros(length(Net)+1,2);
x = cos(theta);
y = sin(theta);
xy(1:length(Net)+1,1) = x(1:length(Net)+1);
xy(1:length(Net)+1,2) = y(1:length(Net)+1);
figure, gplot(Net,xy,'.-');
set(gcf, 'Color', [1 1 1]);
axis('equal');
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
axis off;

format