function [candidates_species, candidates_reactions] = removal_candidate_N (n_s, N, L) 
    %---
    %Candidates for removal
    candidates_species = {};
    candidates_reactions = {};
    n_r = size(L,1) - n_s;
    I = 1:n_s;
    %Adjacecny matrix of the species-reaction graph
    A = diag(diag(L)) - L;
    %All subsets of species of length N
    subsets = nchoosek(I, N);
    %Number of such subsets
    n = size(subsets, 1);
    %---
    for i = 1:n
        %A subset of species for removal
        I_1 = subsets(i,:);
        neighbors_i = [];
        for i1 = 1:numel(I_1)
            i2 = I(i1);
            %Neighbor reactions of species in I_1
            neighbors_i2 = union(find(A(:, i2) ~= 0), find(A(i2, :) ~= 0));
            neighbors_i = union(neighbors_i, neighbors_i2);
        end
        neighbors_i = neighbors_i(:).';
        %Number of these neighbors
        n_i = numel(neighbors_i);
        %All subsets of the neighbor reactions of species in I
        all_subsets_i = arrayfun(@(k) neighbors_i(logical(bitget(k, 1:n_i))), 1:2^n_i-1, 'UniformOutput', false);

        for j = 1:numel(all_subsets_i)
            I_2 = all_subsets_i{j} - n_s;
            if numel(I_2) ~= n_r
                Schur_i = Schur_complement (L, I_1, I_2, n_s);
                if ~isempty(Schur_i) 
                    is_Schur_i_Laplacian = is_Laplacian_matrix (Schur_i, n_s-N);
                    if is_Schur_i_Laplacian
                        I_1 = I_1(:).';
                        I_2 = I_2(:).';
                        candidates_species{end+1} = I_1;
                        candidates_reactions{end+1} = I_2;
                    end
                end   
            end
        end    
    end
   %---
end
