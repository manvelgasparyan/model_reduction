function [Schur] = Schur_complement (L, I_1, I_2, n_s)
   %---
   n_v = size(L, 1);
   %-----
   I = [I_1, I_2 + n_s]; %Indices to delete
   J = setdiff(1:n_v, I); %Indices to remain
   %-----
   L_11 = L(J,J);
   L_12 = L(J,I);
   L_21 = L(I,J);
   L_22 = L(I,I);
   %-----
   d = det(L_22);
   if d ~= 0
      Schur = L_11 - L_12*inv(L_22)*L_21;
   else
       Schur = [];
   end
   %-----
end
