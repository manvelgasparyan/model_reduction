%==========================================================================
%-------The following function corresponds to the complex graph------------ 
%-------of a given biochemical reaction network----------------------------
%==========================================================================
function [Graph] = Complex_Graph(Z,B)
%--------------------------------------------------------------------------
%-------It is a function of two variables Z and B, where Z is the----------
%-------complex composition matrix and B is the incidence matrix-----------
%--------------------------------------------------------------------------
        global s 
        [C,R] = size(B);
        Species = sym('X', [1 s], 'real');
        Species = Species'; %The species of the network
        Complexes = Z'*Species; %The compelxes of the network
        substrates = []; %The tail vertices of the complex graph  
        products = [];   %The corresponding head vertices   
        for i = 1:R 
            for j = 1:C
                for m = 1:C
                    if B(j,i) == -1 && B(m,i) == 1
                       substrates = [substrates, j];
                       products = [products, m];
                    end
                end
            end
        end
        weights = 1:numel(substrates);
        names = string(Complexes);
        Graph = digraph(substrates,products, weights, names);
end
%==========================================================================
