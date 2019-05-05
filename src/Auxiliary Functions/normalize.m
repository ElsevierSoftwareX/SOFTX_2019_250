function [ vector ] = normalize( x )
%NORMALIZE takes a vector x and returns a vector such that the sum of its
%entries is equal to one.
    
    vector = x / sum(x);
end

