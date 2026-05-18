%==========================================================================
%-------The following function determines the set of indinces corresponding
%-------to the species, whose elimination clusters certain reactions-------
%==========================================================================
function [I] = Species_to_Eliminate(Z)
%--------------------------------------------------------------------------
%-------It is a function of one variables Z, which is the complex----------
%-------composition matrix of a given biochemical reaction network---------
%--------------------------------------------------------------------------
global s c
%--------------------------------------------------------------------------
        %Determine the index sets of complexes of the original network
        for i = 1:c
            Index_set{i} = [];
            for j = 1:s
                if Z(j,i) ~= 0
                   Index_set{i} = [Index_set{i},j];  
                end
            end
        end
%--------------------------------------------------------------------------
        %Determine species to eliminate
        Eliminate_species = [];
        no_species = [];
%--------------------------------------------------------------------------        
        for i = 1:c-1
            for j = i+1:c
               Complex_pair = [i, j];
               Shared_species = intersect(Index_set{i},Index_set{j});
               if  numel(Shared_species) ~= 0
                   El_1 = setdiff(Index_set{i},Shared_species);
                   El_2 = setdiff(Index_set{j},Shared_species);
                   El = union(El_1,El_2);
                   Eliminate_species = union(Eliminate_species,El);
               end
               Complex_rest = setdiff(1:c, Complex_pair);
               for n = 1:numel(Complex_rest)
                   m = Complex_rest(n);
                   Remaining{n} = setdiff(Index_set{m}, Eliminate_species);
                   if numel(Remaining{n}) == 0
                      no_species = union(no_species, Index_set{m});
                   end 
               end
       
            end
        end
        Eliminate_species = setdiff(Eliminate_species, no_species);
        I = Eliminate_species; %Species to be eliminated
end
%==========================================================================
