%==========================================================================
%--------The following function is used to determine the steady------------
%--------state of a given biochemical reaction network---------------------
%==========================================================================
function [x_star] = Steady_State(x)
%--------------------------------------------------------------------------
%-------It is a function of one variables x, which is the species'---------
%-------concentration vector-----------------------------------------------
%--------------------------------------------------------------------------
        %The sxc complex composition matrix Z, The cxr incidence matrix B,
        %The vector k of length r of rate constants of the reactions, 
        %The cxr outgoing matrix Delta of B
        global Z B k Delta y_0
%--------------------------------------------------------------------------
        %The diagonal rxr matrix K whose diagonal elements are the elements
        %of k
        K = diag(k);
%--------------------------------------------------------------------------
        %The vector d of length r of rational functions in the expressions 
        %of the reaction rates
        d = Denominator(x);
%--------------------------------------------------------------------------
        %The diagonal rxr matrix D whose diagonal elements are the elements  
        %of d
        D = diag(d);
%--------------------------------------------------------------------------
        %The complex expression vector of length s
        C = exp(Z'*log(x));
%--------------------------------------------------------------------------
       %The vector of reaction rates of length r
       v = -K*D*Delta'*C;
%--------------------------------------------------------------------------       
       %The mxr stoichimetric matrix S
       S = Z*B;       
%--------------------------------------------------------------------------
       %The balance laws 
       dot_x = Z*B*v;
%--------------------------------------------------------------------------
       %The kernel of the transpose matrix of S
       N = null(S');
%--------------------------------------------------------------------------
      %The conservation laws of the network
      Conservation_Laws = N'*x - N'*y_0;
%--------------------------------------------------------------------------
     %We are interested in positive solutions
      Positive = abs(x) - x;
%--------------------------------------------------------------------------
     x_star = [dot_x; Conservation_Laws; Positive];
end
%==========================================================================
