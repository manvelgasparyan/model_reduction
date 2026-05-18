%==========================================================================
%---The following function corresponds to the Laplacian matrix of the------
%---reduced network of a given biochemical reaction network----------------
%==========================================================================
function [Schur] = Laplacian_Reduced(x)
%--------------------------------------------------------------------------
%----It is a function of oe variables x, which is the vector of the--------
%-------------speciec' concentrations of length s--------------------------
%--------------------------------------------------------------------------
      %The (s_eq)x(c_eq) complex composition matrix Z_eq, the (c_eq)x(r_eq)
      %incidence matrix B_eq, the (r_eq)x(c_eq) outward matrix Delta_eq of
      %B_eq, The set I of indices that will be eliminated from the network
        global B_eq Delta_eq I k  
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
        L_eq = B_eq*K*D_eq*Delta_eq';
%--------------------------------------------------------------------------
        %We apply the Rao method to test the reduction
        global Deleted; global Reduced;
%--------------------------------------------------------------------------
        Schur = L_eq(Reduced,Reduced)-L_eq(Reduced,Deleted)*inv(L_eq(Deleted,Deleted))*L_eq(Deleted,Reduced);
%--------------------------------------------------------------------------
 
end
%==========================================================================
