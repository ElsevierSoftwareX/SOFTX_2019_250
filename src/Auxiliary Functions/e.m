function [ vector ] = e( i , n )
%E Creates a canonical vector i of size n.
%   Creates a vector with all elements equal to zero except for the ith
%   entry.
    
    vector = zeros(n,1);
    vector(i,1) = 1;
end

