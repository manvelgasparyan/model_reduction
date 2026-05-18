clear all; close all; clc
%==========================================================================
%==========================================================================
%% PAGE 0: Load model properties
%--------------------------------------------------------------------------
load('data.mat')
%---
assume(species > 0)
%==========================================================================
%==========================================================================
%% PAGE 1: Load model properties
%--------------------------------------------------------------------------
nm = 1:numel(species);
%--------------------------------------------------------------------------
for i = nm
    %----------------------------
    figure;
    %---
    time_1_i = 0:max(tau_1)/n_1:tau_1(i);
    conc_1_interp_i = interp1(time_1, conc_1(:,i), time_1_i);
    %---
    semilogy(time_1_i, conc_1_interp_i, 'r-', 'LineWidth', 1.5); % Blue line, thickness = 2
    %---
    title([char(species(i)) ': Original 1'], 'FontSize', 20, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Number of molecuels (#)', 'FontSize', 18);
    %---
    legend({'Original 1'}, 'FontSize', 18, 'Location', 'best');
    %---
    ax = gca;
    ax.FontSize = 14; % Sets tick labels and axes font size
    ax.XLabel.FontSize = 16;
    ax.YLabel.FontSize = 16;
    %---
    grid on;
    grid minor;
    box on;
    %---
    ax.TickLength = [0.02, 0.02];
    set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]
    %----------------------------
end
%--------------------------------------------------------------------------
for i = nm
    %----------------------------
    figure;
    %---
    time_2_i = 0:max(tau_2)/n_2:tau_2(i);
    conc_2_interp_i = interp1(time_2, conc_2(:,i), time_2_i);
    %---
    semilogy(time_2_i, conc_2_interp_i, 'b-', 'LineWidth', 1.5); % Blue line, thickness = 2
    %---
    title([char(species(i)) ': Original 2'], 'FontSize', 20, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Number of molecuels (#)', 'FontSize', 18);
    %---
    legend({'Original 2'}, 'FontSize', 16, 'Location', 'best');
    %---
    ax = gca;
    ax.FontSize = 14; % Sets tick labels and axes font size
    ax.XLabel.FontSize = 16;
    ax.YLabel.FontSize = 16;
    %---
    grid on;
    grid minor;
    box on;
    %---
    ax.TickLength = [0.02, 0.02];
    set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]
    hold off;
    %----------------------------
end
%--------------------------------------------------------------------------
%{
    %---
    figure;
    %Plot the graph
    [n_s,n_r] = size(stoich);
    pG = plot(sr_graph, 'Layout', 'force');
    pG.EdgeFontSize = 7; 
    %Remove the initial labels
    pG.NodeLabel = {};
    %Set label positions manually
    for i = 1:numel(pG.XData)
       text(pG.XData(i), pG.YData(i), sr_graph.Nodes.Name{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    %Customize the appearance of species
    highlight(pG, 1:n_s, 'NodeColor', [0.0 0.6 1.0], 'Marker', 'o', 'MarkerSize', 20);
    %Customize the appearance of reactions
    highlight(pG, n_s+1:n_s+n_r, 'NodeColor', [1 0.8 0.0], 'Marker', 's', 'MarkerSize', 20);
    %Customize the appearance of edges
    set(pG, 'EdgeColor', [0, 0, 0], 'LineWidth', 1, 'ArrowSize', 5, 'EdgeAlpha', 1);
    %---
end
%}
%==========================================================================
%==========================================================================
