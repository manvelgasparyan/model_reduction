clear all; close all; clc
%==========================================================================
%==========================================================================
%% PAGE 0: Load model properties
%--------------------------------------------------------------------------
load('data.mat')
%---
assume(species > 0)
assume(p > 0)
assume(b > 0)
%==========================================================================
%==========================================================================
%% PAGE 1: Laplacian model
%--------------------------------------------------------------------------
weights = sr_graph_weights (stoich, rates, species);
%----------
L = Laplacian_matrix(weights, sr_graph);
%----------
[d, Lambda, v] = Laplacian_model (species, stoich, L);
%----------
[v_sol, ds_dt_sol] = Laplacian_ODE(L, v, d, rates, species);
%==========================================================================
%==========================================================================
%% PAGE 2: Save the output
%--------------------------------------------------------------------------
save('data.mat', 'weights', '-append');
save('data.mat', 'L', '-append');
save('data.mat', 'd', '-append');
save('data.mat', 'Lambda', '-append');
save('data.mat', 'v', '-append');
save('data.mat', 'v_sol', '-append');
save('data.mat', 'ds_dt_sol', '-append');
%==========================================================================
%% PAGE 3: Main function
%--------------------------------------------------------------------------
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
               weights(weight_count) = rates(j)/species(i);
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
   weights = weights(1:weight_count);
end
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
function [d, Lambda, v] = Laplacian_model (species, stoich, L)
%---
    d = - sum(stoich,1);
    %---
    n_r = size(stoich,2);
    %---
    v = sym('v', [n_r,1]);
    assume(v > 0)
    %---
    s = species';
    %---
    Lambda = - L*[s; v];
%---
end
%--------------------------------------------------------------------------
function [v_sol, ds_dt_sol] = Laplacian_ODE(L, v, d, rates, species)
    %---
    n_r = length(v);
    n_s = length(species);
    %---
    L_11 = L(1:n_s,1:n_s);
    L_12 = L(1:n_s,1+n_s:n_r+n_s);
    L_21 = L(1+n_s:n_r+n_s,1:n_s);
    L_22 = L(1+n_s:n_r+n_s,1+n_s:n_r+n_s);
    %---
    Gamma_d = diag(d);
    %---
    v_sol = -(pinv(L_22+Gamma_d)*L_21)*species';
    %---
    for j = 1:n_r
        if v_sol(j) ==0
           v_sol(j) = rates(j);
        end
    end
    %---
    ds_dt_sol = -L_11*species' - L_12*v_sol;
    %---
end
%==========================================================================
