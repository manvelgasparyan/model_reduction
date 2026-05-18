function [G] = sr_graph_L (n_s, L, I_1, I_2)
   %---
   n_r = size(L, 1) - n_s;
   %---
   J_1 = setdiff(1:n_s, I_1);
   J_2 = setdiff(1:n_r, I_2);
   %---
   %Create a directed bipartite graph
   G = digraph();
   %---
   %Add vertices for species
   for i1 = 1:numel(J_1)
       j1 = J_1(i1);
       G = addnode(G, ['S_{' num2str(j1), '}']);
   end
   %Add vertices for reactions
   for i2 = 1:numel(J_2)
       j2 = J_2(i2);
       G = addnode(G, ['R_{' num2str(j2), '}']);
   end
   %---
   Schur = Schur_complement (L, I_1, I_2, n_s);
   %---
   n_v = size(Schur, 1);
   %---
   edges = [];
   for i = 1:n_v
       for j = 1:n_v
           if i ~= j && Schur(i,j) ~= 0
               edge_i = [j, i];
               edges = [edges; edge_i];
           end
       end
   end
   %---
   %Add the edges to the graph
   for i = 1:size(edges, 1)
       source_node = edges(i, 1);
       target_node = edges(i, 2);
       G = addedge(G, source_node, target_node); 
   end
end
