%==========================================================================
%---The following function corresponds to the Rao test model of the given--
%---biochemical reaction network-------------------------------------------
%==========================================================================
function [dot_x] = Rao_Model(~, x)
%--------------------------------------------------------------------------
%----It is a function of one variables x, which is the vector of the-------
%-------------species' concentrations of length s--------------------------
%--------------------------------------------------------------------------
      %The (s_eq)x(c_eq) complex composition matrix Z_eq, the (c_eq)x(r_eq)
      %incidence matrix B_eq, the (r_eq)x(c_eq) outward matrix Delta_eq of
      %B_eq, The set I of indices that will be eliminated from the network
       global Z_tilde B_tilde Delta_tilde I k
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
        L_eq = B_tilde*K*D_eq*Delta_tilde';
%--------------------------------------------------------------------------
        %We apply the Rao method to test the reduction
        global del red
%--------------------------------------------------------------------------
        %The complex expression vector of length s
        C_tilde = exp(Z_tilde(:,red)'*log(x)); 
%--------------------------------------------------------------------------
       Schur = L_eq(red,red)-L_eq(red,del)*inv(L_eq(del,del))*L_eq(del,red);
%--------------------------------------------------------------------------
        dot_x = -Z_tilde(:,red)*Schur*C_tilde;
end
%==========================================================================
