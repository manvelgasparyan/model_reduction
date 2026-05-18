function [L] = Laplacian_matrix(weights, sr_graph)
%---
   Gamma_weights = diag(weights);
   %---
   B = full(incidence(sr_graph));
   %---
   Delta = min(B,0);
   %---
   L = B*Gamma_weights*Delta';
%---
end
