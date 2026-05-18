function [weights] = sr_graph_weights (stoich, rates, species)
%----------------
   [n_s,n_r] = size(stoich); 
   %Maximum number of edges
   n_e = n_s * n_r; 
   %Initialize edges matrix  
   weights = sym(zeros(n_e, 2)); 
   %Counter for the number of edges added
   weight_count = 0;
   %Iterate over the matrix S
   for i = 1:n_s
       for j = 1:n_r
           S_ij = stoich(i,j);
           %Check if S_ij is negative
           if S_ij < 0
               weight_count = weight_count + 1;
               weights(weight_count) = diff(rates(j),species(i));
           end
       end
   end
   %---
   for i = 1:n_s
       for j = 1:n_r
           S_ij = stoich(i,j);
           %Check if S_ij is positive
           if S_ij > 0
              weight_count = weight_count + 1;
              weights(weight_count) = sym(S_ij);
           end
       end
   end
   %Remove excess rows if not all edges are used
   weights = weights(1:weight_count)';
end
