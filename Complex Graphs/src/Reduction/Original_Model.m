%==========================================================================
%-------The following function corresponds to the original mathematical---- 
%-------model of a given biochemical reaction network----------------------
%==========================================================================
function [dot_x] = Original_Model(~, x)
%--------------------------------------------------------------------------
%-------It is a function of two variables ~ and x, where x is the vector---
%-------of species' concentrations and ~ is the time variable--------------
%--------------------------------------------------------------------------
        %The sxc complex composition matrix Z, The cxr incidence matrix B,
        %The vector k of length r of rate constants of the reactions, The 
        %cxr outgoing matrix Delta of B
        global Z B k Delta 
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
        %The cxc Laplacian matrix 
        L = B*K*D*Delta';
%--------------------------------------------------------------------------
        %The complex expression vector of length s
        C = exp(Z'*log(x));
%--------------------------------------------------------------------------
        %The balance laws representing the tare of change of the 
        %concentration functions of s species
        dot_x = -Z*L*C;
end
%==========================================================================
