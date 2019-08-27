function w_max = optQuasiConvex(A,functionToGetTransition,G)
%OPTQUASICONVEX Returns the optimal w that guarantees the fastest convergence for
%any algorithm with quasiconvex property on the optimization of T that might depend on G.
%   Solves the optimization problem:
%       min    max{ lambda( T ) }
%   Given the format of the problem, the solution is found through
%   bissetion.

lower = 0;
upper = 2;
epsilon = 1E-5;

while (upper-lower) >= epsilon
    w = (lower+upper)/2;
    
    Tw = functionToGetTransition(A,G,w);
    dTw = functionToGetTransition(A,G,w+epsilon);
    
    if spectralRadius(dTw) < spectralRadius(Tw)
        lower = w+epsilon;
    else
        upper = w;
    end
end
w_max = (lower+upper)/2;
end

