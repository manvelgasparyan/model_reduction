%==========================================================================
%--The following function determines the suitable conservation laws that---
%-----will be used to elliminate the species with the set of indices I-----
%==========================================================================
function [Elimination] = Elimination_of_species(x)
%----It is a function of one variables x, which is the vector of the-------
%-------------species' concentration functions of length s-----------------
%--------------------------------------------------------------------------
        %The sxc complex composition matrix Z, The cxr incidence matrix B,
        %The vector k of length r of proportionality constants of the
        %reactions, The set I of indices that will be eliminated from 
        %the network
        global Z B y_0 r I
%--------------------------------------------------------------------------
        %The sxr stoichiometric matrix S of the network 
        S = Z*B;
%--------------------------------------------------------------------------
        %For eversy species that we want to eliminate from the network we
        %find the suitable conservation law
        for j = 1:numel(I)
            p = I(j);
            I1 = setdiff(I,p);
            n1 = numel(I1);
            S1 = S;
            S1(I1,:) = zeros(n1,r);
            N = null(S1','r');
            [~,cl] = size(N); 
            for i = 1:cl
                if N(p,i) ~= 0
                    M = N(:,i);
            break
                end
            end
            CL(:,j) = M;
            C = M'*y_0;
            CQ(j) = C;
            a = M(p);
            M(p) = 0;
            Elimination(j) = (-M'*x+C)/a;
        end
end
%==========================================================================
