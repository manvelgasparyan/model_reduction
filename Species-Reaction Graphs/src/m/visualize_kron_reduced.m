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
%% PAGE 1: Visualize
%--------------------------------------------------------------------------
for i = 1:numel(species_red)
    %----------------------------
    figure;
    %---
    time_1_i = 0:max(tau_1_red)/n_1:tau_1_red(i);
    conc_1_interp_i = interp1(time_red_1, conc_red_1(:,i), time_1_i);
    %---
    plot(time_1_i, conc_1_interp_i, 'r--', 'LineWidth', 1.5); % Blue line, thickness = 2
    %---
    title([char(species(i)) ': Reduced 1'], 'FontSize', 20, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Number of molecuels (#)', 'FontSize', 18);
    %---
    legend({'Reduced 1'}, 'FontSize', 18, 'Location', 'best');
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
for i = 1:numel(species_red)
    %----------------------------
    figure;
    %---
    time_2_i = 0:max(tau_2_red)/n_2:tau_2_red(i);
    conc_2_interp_i = interp1(time_red_2, conc_red_2(:,i), time_2_i);
    %---
    plot(time_2_i, conc_2_interp_i, 'b--', 'LineWidth', 1.5); % Blue line, thickness = 2
    %---
    title([char(species(i)) ': Reduced 2'], 'FontSize', 20, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Number of molecuels (#)', 'FontSize', 18);
    %---
    legend({'Reduced 2'}, 'FontSize', 18, 'Location', 'best');
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

%==========================================================================
%==========================================================================
