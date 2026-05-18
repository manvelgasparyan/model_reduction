function [G] = sr_graph (stoich)
   %---
   %The number of species and reactions
   [n_s,n_r] = size(stoich);
   %---
   %Create a directed bipartite graph
   G = digraph();
   %---
   %Add vertices for species
   for i = 1:n_s
       G = addnode(G, ['S_{' num2str(i), '}']);
   end
   %Add vertices for reactions
   for j = 1:n_r
       G = addnode(G, ['R_{' num2str(j), '}']);
   end
   %---
   %Maximum number of edges
   n_e = n_s * n_r; 
   %Initialize edges matrix  
   edges = zeros(n_e, 2); 
   %Counter for the number of edges added
   edge_count = 0;
   %Iterate over the matrix S
   for i = 1:n_s
       for j = 1:n_r
           S_ij = stoich(i,j);
           %Check if S_ij is negative
           if S_ij < 0
               edge_count = edge_count + 1;
               edges(edge_count, :) = [i, n_s+j];
           end
       end
   end
   %---
   for i = 1:n_s
       for j = 1:n_r
           S_ij = stoich(i,j);
           %Check if S_ij is positive
           if S_ij > 0
              edge_count = edge_count + 1;
              edges(edge_count, :) = [n_s+j, i];
           end
       end
   end
   %Remove excess rows if not all edges are used
   edges = edges(1:edge_count, :);
   %---
   %Add the edges to the graph
   for i = 1:size(edges, 1)
       source_node = edges(i, 1);
       target_node = edges(i, 2);
       G = addedge(G, source_node, target_node); 
   end
   %---
end
