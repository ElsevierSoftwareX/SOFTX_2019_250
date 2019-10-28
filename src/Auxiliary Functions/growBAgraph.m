function A = growBAgraph(N,S,m)
% Generates a scale-free directed adjacency matrix using Barabasi and
% Albert algorithm from an existing network.
% Input: N - number of nodes in the network, S - seed network, m = average degree.
% Example: A = BAgraph(300,previousA,10);
% Ref: Methods for generating complex networks with selected structural properties for simulations, Pretterjohn et al, Frontiers Comp Neurosci

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

hwait = waitbar(0,'Please wait. Growing a directed scale-free adjacency matrix');
A = zeros(N);

% Compute size of the seed matrix S
mo = length(S);

% Copy the seed network to the new adjacency matrix
A(1:mo,1:mo) = S;

% Update current number of links
E = sum(sum(S));

% Second add remaining nodes with a preferential attachment bias - rich get
% richer
for i=mo+1:N
	waitbar(i/N,hwait,sprintf('Please wait. Generating directed scale-free adjacency matrix\n%.1f %%',(i/N)*100));
	curr_deg =0;
	while(curr_deg<m)
		sample = setdiff(1:N,[i find(A(i,:))]);
		j = datasample(sample,1);
		b = sum(A(j,:))/E;
		r = rand(1);
		if(b>r)
			r = rand(1);
			if(b>r)
			A(i,j) = 1;
			A(j,i) = 1;
			E = E +2;
			else
			A(i,j) = 1;
			E = E +1;
			end
		else
			no_connection = 1;
			while(no_connection)
				sample = setdiff(1:N,[i find(A(i,:))]);
				h = datasample(sample,1);
				b = sum(A(h,:))/E;
				r = rand(1);
				if(b>r)
					r = rand(1);
					if(b>r)
					A(h,i) = 1;
					A(i,h) = 1;
					E = E +2;
					else
					A(i,h) = 1;
					E = E+1;
					end
					no_connection = 0;
				end
			end
		end
		curr_deg = sum(A(i,:));
	end
end
delete(hwait);

