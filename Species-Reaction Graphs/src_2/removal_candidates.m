function [candidates_species, candidates_reactions] = removal_candidates (L, n_s)
%---
candidates_species = {};
candidates_reactions = {};
%---
for N = 1:n_s-2
    [candidates_species_N, candidates_reactions_N] = removal_candidate_N (n_s, N, L);
    candidates_species = [candidates_species, candidates_species_N];
    candidates_reactions = [candidates_reactions, candidates_reactions_N];
%---
end
