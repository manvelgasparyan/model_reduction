%==========================================================================
%--------The following function corresponds to the rational functions------
%------------in the expressions of the reaction rates in the---------------
%--------------------equivalent reaction network.--------------------------
%==========================================================================
function [d_eq] = Denominator_Equivalent(x)
%--------------------------------------------------------------------------
%----It is a function of one variables x, which is the vector of the-------
%-------------species' concentrations of length s--------------------------
%--------------------------------------------------------------------------
        %The sxc complex composition matrix Z, The cxr incidence matrix B,
        %The vector k of length r of proportionality constants of the
        %reactions, The cxr outgoing matrix Delta, The set I of indices 
        %that will be eliminated from the network
        global Z B I r
%--------------------------------------------------------------------------
        %The sxr stoichiometric matrix S of the network 
        S = Z*B;
%--------------------------------------------------------------------------
        d = Denominator(x);
        for i =1:numel(I)
            ind = I(i);
            for j =1:r
                if S(ind,j)==-1 
                   d(j) = d(j)*x(ind);
                end
            end
        end
        d_eq = d;
end
%==========================================================================
