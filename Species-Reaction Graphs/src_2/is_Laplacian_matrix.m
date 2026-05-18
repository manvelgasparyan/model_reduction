function [result] = is_Laplacian_matrix(L, n_s)
   %---
   n_r = size(L,1) - n_s;
   %---
   %Species indices
   I_1 = 1:n_s;
   %Reaction indices
   I_2 = n_s+1:n_s+n_r;
   %---
   %Species block
   L_11 = L(I_1,I_1);
   %Reactions-Species block
   L_12 = L(I_1, I_2);
   %Species-Reactions block
   L_21 = L(I_2,I_1);
   %Reactions block
   L_22 = L(I_2, I_2);
   %---
   is_11_diagonal = isequal(L_11, diag(diag(L_11)));
   is_22_diagonal = isequal(L_22, diag(diag(L_22)));
   is_22_constant = min(min(isAlways(~has(L_22, symvar(L_22)))));
   is_12_constant = min(min(isAlways(~has(L_12, symvar(L_12)))));
   is_not_bidirectional = is_not_reverse_edges(L, n_s);
   %---
   if is_11_diagonal && is_22_diagonal && is_22_constant && is_12_constant && is_not_bidirectional
     result = 1;
     %fprintf("The given matrix is a Laplacian of an SR-graph")
   else
     result = 0;
     %fprintf("The given matrix is not a Laplacian of an SR-graph")
   end
   %-------------------------------------------
   %-------------------------------------------
   function [result] = is_not_reverse_edges(L, n_s)
       %----------------
       value = 0;
       n_v = size(L,1);
       for i = 1:n_s
           for j = n_s+1:n_v
               value = value + abs(L(i,j)*L(j,i));
           end
       end
       %----------------
       if value == 0
           result = 1;
       else
           result = 0;
       end
       %----------------
   
   end
   %-------------------------------------------
end
