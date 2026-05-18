%==========================================================================
%---The following function corresponds to the reduced model of a given-----
%---biochemical reaction network.
%==========================================================================
function [dot_x] = Reduced_Model(~, x)
%--------------------------------------------------------------------------
%----It is a function of one variables x, which is the vector of the-------
%-------------species' concentrations of length s--------------------------
%--------------------------------------------------------------------------
      %The (s_eq)x(c_eq) complex composition matrix Z_eq, The set I of 
       %indices that will be eliminated from the network
       global Z_eq I 
%--------------------------------------------------------------------------
        Elimination = Elimination_of_species(x);
        for j = 1:numel(I)
            p = I(j);
            x(p) = Elimination(j);
        end        
%--------------------------------------------------------------------------
        %We apply the Rao method to test the reduction
        global Reduced;
%--------------------------------------------------------------------------
        %The complex expression vector of length s
        C_eq = exp(Z_eq(:,Reduced)'*log(x)); 
%--------------------------------------------------------------------------
       Schur = Laplacian_Reduced(x);
%--------------------------------------------------------------------------
        dot_x = -Z_eq(:,Reduced)*Schur*C_eq;
end
%==========================================================================
