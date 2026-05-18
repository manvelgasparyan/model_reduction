%==========================================================================
%-------The following function corresponds to the equivalent model---------
%---------------of a given biochemical reaction network.-------------------
%==========================================================================
function [dot_x] = Equivalent_Model(~, x)
%--------------------------------------------------------------------------
%----It is a function of one variables x, which is the vector of the-------
%-------------species' concentrations of length s--------------------------
%--------------------------------------------------------------------------
      %The (s_eq)x(c_eq) complex composition matrix Z_eq, the (c_eq)x(r_eq)
      %incidence matrix B_eq, the (r_eq)x(c_eq) outgoing matrix Delta_eq
      %The set I of indices that will be eliminated from the network
       global Z_eq B_eq Delta_eq k I
%--------------------------------------------------------------------------
        Elimination = Elimination_of_species(x);
        for j = 1:numel(I)
            p = I(j);
            x(p) = Elimination(j);
        end        
%--------------------------------------------------------------------------
        %The diagonal rxr matrix K whose diagonal elements are the elements
        %of k
        K = diag(k);
%--------------------------------------------------------------------------
        %The vector d of length r of rational functions in the expressions 
        %of the reaction rates
        d_eq = Denominator_Equivalent(x);
%--------------------------------------------------------------------------
        %The diagonal rxr matrix D_eq whose diagonal elements are the  
        %elements of d_eq
        D_eq = diag(d_eq);        
%--------------------------------------------------------------------------
        %The Laplacian matrix of the equivalent network
        L_eq = B_eq*K*D_eq*Delta_eq';
%--------------------------------------------------------------------------
        %The complex expression vector of length s
        C_eq = exp(Z_eq'*log(x));        
%--------------------------------------------------------------------------
        %The balance laws of the equivalent network
        dot_x = -Z_eq*L_eq*C_eq;
end
%==========================================================================
