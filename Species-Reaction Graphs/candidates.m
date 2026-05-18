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
%% PAGE 1: Candidates
%--------------------------------------------------------------------------
n_s = length(species);
n_r = length(rates);
%---
I_s = [9, 11, 12, 13, 14];
I_r = [5:9, 11:15, 17, 18];
%---
sr_graph_red = remove_species_reactions (sr_graph, I_s, I_r, species)
%---
species_nodes = 1:n_s;
reactions_nodes = 1:n_r;
result = is_sr_graph (sr_graph, species_nodes, reactions_nodes)
%==========================================================================
%==========================================================================
%% PAGE 2: Main functions
%--------------------------------------------------------------------------
function [sr_graph_red] = remove_species_reactions (sr_graph, I_s, I_r, species)
%---
    N = numnodes(sr_graph);
    %---
    n_s = length(species);
    %---
    G = sr_graph;
    %---
    I = [I_s,I_r + n_s];
    %---
    for i1 = 1:length(I)
        i2 = I(i1);
        for j1 = 1:N
            for j2 = j1:N
                e1 = findedge(sr_graph, j1, i2) > 0;
                e2 = findedge(sr_graph, i2, j2) > 0;
                if e1 == 1 && e2 ==1
                   G = addedge(G, j1, j2); 
                end
            end
        end
    end
    %---
    G = rmnode(G, I);
    %---
    sr_graph_red = G;
%---
end
%--------------------------------------------------------------------------
function [result] = is_sr_graph (di_graph, species_nodes, reactions_nodes)
%---
    A = adjacency(di_graph); 
%---
    I_1 = species_nodes;
    I_2 = reactions_nodes + length(species_nodes);
%---
    subA_1 = A(I_1, I_1);          
    result_1 = any(subA_1(:)); %true if no edges from I_1 to itself
%---
    subA_2 = A(I_2, I_2);          
    result_2 = any(subA_2(:)); %true if no edges from I_2 to itself
%---
    d = indegree(di_graph) + outdegree(di_graph);
    result_3 = all(d > 0);   %true → every node is connected, and false → there is at least one isolated node.
%---
result = hasEdgesWithin;
%---
end
%==========================================================================
%==========================================================================
