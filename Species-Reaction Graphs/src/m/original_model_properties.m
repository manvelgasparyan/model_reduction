clear all; close all; clc
%==========================================================================
%==========================================================================
%% PAGE 0: Model properties
%--------------------------------------------------------------------------
syms Anchor PP1 MR MR1 MR2 C5 C6 IR IR1 IR2 C1 C2 C3 C4
species = [Anchor, PP1, MR, MR1, MR2, C5, C6, IR, IR1, IR2, C1, C2, C3, C4];
assume(species > 0)
%==========================================================================
%==========================================================================
%% PAGE 1: Original model
%--------------------------------------------------------------------------
init_cond_1 = [167.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
%---
%init_cond_2 = [0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 165.0, 0.0];
init_cond_2 = [0.0, 2.0, 167.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
%--------------------------------------------------------------------------
%Parameters
theta = [0.0002; 0.008; 0.0008; 0.018; 1.0; 2e-05; 2e-05; 0.14815; 24.0; 
         6.0; 0.14815; 24.0; 6.0; 0.14815; 24.0; 6.0; 0.14815; 24.0; 6.0;
         0.29167; 1.4; 0.35; 0.29167; 1.4; 0.35];
%--------------------------------------------------------------------------
syms k1f k1b k2 k3f k3b k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 
syms k15 k16 k17 k18 k19 k20 k21 k22 k23 AMPAR PKA 
%--------------------------------------------------------------------------
rates = [k1f*Anchor*IR2, k1b*MR2, k2*MR, k3f*AMPAR, k3b*IR, k4*IR2, ...
         k5*IR1, k6*PKA*MR, k7*C1, k8*C1, k9*PKA*IR, k10*C2, k11*C2, ...
         k12*PKA*MR1, k13*C3, k14*C3, k15*PKA*IR1, k16*C4, k17*C4, ...
         k18*MR1*PP1, k19*C5, k20*C5, k21*MR2*PP1, k22*C6, k23*C6];
%--------------------------------------------------------------------------
p = [k1f; k1b; k2; k3f; k3b; k4; k5; k6; k7; k8; k9; k10; k11; k12; k13;
     k14; k15; k16; k17; k18; k19; k20; k21; k22; k23];
%---
assume(p > 0)
%--------------------------------------------------------------------------
b = [AMPAR; PKA]; 
%---
assume(b > 0)
%--------------------------------------------------------------------------
%Stoichiometric matrix
stoich = [-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 1, 1;
 0, 0, -1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1;
 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1;
 0, 0, 1, 1, -1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0;
 -1, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0];
%==========================================================================
%==========================================================================
%% PAGE 3: Integration
%--------------------------------------------------------------------------
tau_1 = [2*10^4; 2*10^4; 2*10^4; 2*10^4; 1.7*10^4; 1.2*10^4; 1.2*10^4; 7; 300; 1.5*10^4; 1.2*10^4; 10; 1.2*10^4; 200];
%---
tau_2 = [2.6*10^5; 0.6*10^5; 230; 150; 4.5*10^4; 150; 150; 10; 200; 2.2*10^5; 300; 7; 0.3; 150];
%--------------------------------------------------------------------------
tau_1 = 1.1*tau_1;
tau_2 = 1.1*tau_2;
%--------------------------------------------------------------------------
m = 7;
%--------------------------------------------------------------------------
n_1 = 10^m;
%---
n_2 = 10^m;
%--------------------------------------------------------------------------
ds_dt = @(t, s)original_ode_model (t, s, theta, stoich);
%--------------------------------------------------------------------------
tic
[time_1, conc_1] = ODE_integration (ds_dt, init_cond_1, tau_1, n_1);
elapsed_time = toc;
disp(['Elapsed time (original 1): ', num2str(elapsed_time), ' seconds']);
%---
tic
[time_2, conc_2] = ODE_integration (ds_dt, init_cond_2, tau_2, n_2);
elapsed_time = toc;
disp(['Elapsed time (original 2): ', num2str(elapsed_time), ' seconds']);
%==========================================================================
%==========================================================================
%% PAGE 4: Species-Reaction graph
%--------------------------------------------------------------------------
sr_graph = species_reaction_graph (stoich, species);
%==========================================================================
%==========================================================================
%% PAGE 5: Save the output
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
save('data.mat', 'species');
save('data.mat', 'init_cond_1', '-append');
save('data.mat', 'init_cond_2', '-append');
save('data.mat', 'theta', '-append');
save('data.mat', 'p', '-append');
save('data.mat', 'b', '-append');
save('data.mat', 'rates', '-append');
save('data.mat', 'stoich', '-append');
save('data.mat', 'tau_1', '-append');
save('data.mat', 'tau_2', '-append');
save('data.mat', 'n_1', '-append');
save('data.mat', 'n_2', '-append');
save('data.mat', 'time_1', '-append');
save('data.mat', 'time_2', '-append');
save('data.mat', 'conc_1', '-append');
save('data.mat', 'conc_2', '-append');
save('data.mat', 'sr_graph', '-append');
%==========================================================================
%==========================================================================
%% PAGE 6: Functions
%--------------------------------------------------------------------------
function [time, conc] = ODE_integration (ds_dt, init_cond, tau, n)
    %---
    T = max(tau);
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);
    %---
    [time, conc] = ode15s(ds_dt, 0:T/n:T, init_cond, options);
    %---
end
%--------------------------------------------------------------------------
function [ds_dt] = original_ode_model (~, s, theta, stoich)
   %---
   AMPAR = 9;
   PKA = 1;   
   %---
   k1f = theta(1); k1b = theta(2); k2 = theta(3);  k3f = theta(4);
   k3b = theta(5); k4 = theta(6);  k5 = theta(7);  k6 = theta(8);
   k7 = theta(9);  k8 = theta(10); k9 = theta(11); k10 = theta(12);
   k11 = theta(13); k12 = theta(14); k13 = theta(15); k14 = theta(16);
   k15 = theta(17); k16 = theta(18); k17 = theta(19); k18 = theta(20);
   k19 = theta(21); k20 = theta(22); k21 = theta(23); k22 = theta(24);
   k23 = theta(25);   
   %---
   Anchor = s(1);       PP1 = s(2);         MR = s(3); 
   MR1 = s(4);          MR2 = s(5);         C5 = s(6); 
   C6 = s(7);           IR = s(8);          IR1 = s(9); 
   IR2 = s(10);         C1 = s(11);         C2 = s(12);
   C3 = s(13);          C4 = s(14);
   %---
   rates = [k1f*Anchor*IR2, k1b*MR2, k2*MR, k3f*AMPAR, k3b*IR, k4*IR2, ...
         k5*IR1, k6*PKA*MR, k7*C1, k8*C1, k9*PKA*IR, k10*C2, k11*C2, ...
         k12*PKA*MR1, k13*C3, k14*C3, k15*PKA*IR1, k16*C4, k17*C4, ...
         k18*MR1*PP1, k19*C5, k20*C5, k21*MR2*PP1, k22*C6, k23*C6];
   %---
   ds_dt = stoich*rates';   
   %---
end
%--------------------------------------------------------------------------
function [G] = species_reaction_graph (stoich, species)
   %---
   %The number of species and reactions
   [n_s,n_r] = size(stoich);
   %---
   %Create a directed bipartite graph
   G = digraph();
   %---
   %Add vertices for species
   for i = 1:n_s
       G = addnode(G, char(species(i)));
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
%==========================================================================
%==========================================================================
