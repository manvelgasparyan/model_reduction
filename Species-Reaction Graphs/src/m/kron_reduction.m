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
%% PAGE 1: Kron reduction
%--------------------------------------------------------------------------
I_s = [9, 11, 12, 13, 14];  %species to delete
%---
I_r = [5:9, 11:15, 17, 18]; %reactions to delete
%--------------------------------------------------------------------------
tic
[Schur, d_red, v_red, species_red, Lambda_red] = Kron_reduction (species, v, L, d, I_s, I_r);
elapsed_time = toc;
disp(['Elapsed time (Kron reduction): ', num2str(elapsed_time), ' seconds']);
%--------------------------------------------------------------------------
[v_red, ds_dt_red] = Kron_reduced_Laplacian_ODE(Schur, v_red, d_red, species_red, rates);
%--------------------------------------------------------------------------
ds_dt_red = subs(ds_dt_red, p, theta);
%ds_dt_red = subs(ds_dt_red, b(2), 1);
ds_dt_red = simplify(ds_dt_red);
%--------------------------------------------------------------------------
v_red_0 = subs(v_red, p, theta);
v_red_0 = subs(v_red_0, b(2), 1);
%v_red = subs(v_red, b(2), 1);
v_red_0 = simplify(v_red_0);
%==========================================================================
%==========================================================================
%% PAGE 2: Save the output
%--------------------------------------------------------------------------
save('data.mat', 'Schur', '-append');
save('data.mat', 'd_red', '-append');
save('data.mat', 'v_red', '-append');
save('data.mat', 'species_red', '-append');
save('data.mat', 'Lambda_red', '-append');
save('data.mat', 'v_red', '-append');
save('data.mat', 'ds_dt_red', '-append');
%==========================================================================
%% PAGE 3: Main function
%--------------------------------------------------------------------------
function [Schur,d_red, v_red, species_red, Lambda_red] = Kron_reduction (species, v, L, d, I_s, I_r)
   %-----
   n_s = length(species);
   n_r = length(v);
   %-----
   J_s = setdiff(1:n_s, I_s); %Species to remain
   %---
   J_r = setdiff(1:n_r, I_r); %Reactions to remain
   %-----
   I_sr = [I_s, I_r + n_s]; %Indices to delete
   J_sr = setdiff(1:n_s+n_r, I_sr); %Indices to remain
   %-----
   L_11 = L(J_sr,J_sr);
   L_12 = L(J_sr,I_sr);
   L_21 = L(I_sr,J_sr);
   L_22 = L(I_sr,I_sr);
   %-----
   Schur = L_11 - L_12*pinv(L_22)*L_21;
   %-----
   d_red= d(J_r);
   %-----
   species_red = species(J_s);
   %-----
   v_red = v(J_r);
   %-----
   Lambda_red = - Schur*[species_red'; v_red];
   %-----
end
%--------------------------------------------------------------------------
function [v_sol, ds_dt_sol] = Kron_reduced_Laplacian_ODE(L, v, d, species, rates)
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
%--------------------------------------------------------------------------
