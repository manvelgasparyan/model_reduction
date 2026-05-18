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
%---
n_s = length(species);
I_s = [9, 11, 12, 13, 14];  %species to delete
J_s = setdiff(1:n_s, I_s); %Species to remain
%---
conc_1_contracted = conc_1(:,J_s);
conc_2_contracted = conc_2(:,J_s);
%==========================================================================
%==========================================================================
%% PAGE 1: Visualize
%--------------------------------------------------------------------------
for i = 1:numel(species_red)
    %----------------------------
    fig = figure;
    fig.Units = 'pixels';
    % fig.Position = [100, 100, 600, 200];  % This controls outer figure size

    % Now fix the axes position manually:
    % ax = axes('Parent', fig);
    % ax.Units = 'normalized';
    % ax.Position = [0.13, 0.15, 0.775, 0.75];  % [left, bottom, width, height] as a fraction of fig
    
    time_1_i = 0:max(tau_1)/n_1:tau_1(i);
    %time_1_i = time_1;
    %---
    conc_1_interp_i = interp1(time_1, conc_1_contracted(:,i), time_1_i);
    conc_1_red_interp_i = interp1(time_1, conc_red_1(:,i), time_1_i);
    %---
    conc_2_interp_i = interp1(time_1, conc_2_contracted(:,i), time_1_i);
    conc_2_red_interp_i = interp1(time_1, conc_red_2(:,i), time_1_i);
    %---
    semilogy(time_1_i, conc_1_interp_i, '-', 'Color', [0, 0, 0], 'LineWidth', 3); % Org 1 
    hold on
    semilogy(time_1_i, conc_1_red_interp_i, '-.', 'Color',  [0.902, 0.427, 0.055], 'LineWidth', 4); % Red. 1 
    %---
    hold on
    semilogy(time_1_i, conc_2_interp_i, '-', 'Color', [0 0.4470 0.7410],  'LineWidth', 3); % Org 2
    hold on
    semilogy(time_1_i, conc_2_red_interp_i, '-.', 'Color', [1 0 0], 'LineWidth', 4); % Red. 2
    %---
    hold off
    %---
    title(char(species_red(i)), 'FontSize', 25, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 20);
    %     ylabel('Number of molecuels (log scale, unit #)', 'FontSize', 20);
    ylabel(sprintf('Number of molecules \n (log scale, unit #)'), 'FontSize', 20);

    %---
    legend({'Original 1', 'Reduced 1', 'Original 2', 'Reduced 2'}, 'FontSize', 20, 'Location', 'best');
    %---
    ax = gca;
    ax.FontSize = 18; % Sets tick labels and axes font size
    ax.XLabel.FontSize = 18;
    ax.YLabel.FontSize = 18;
    %---
    %grid on;
    %grid minor;
    box on;
    %---
    ax.TickLength = [0.02, 0.02];
    set(gcf, 'Position', [100, 100, 800, 450]); % [left, bottom, width, height]
    %ylim([0 10^3]);  % This is fixed y-lim, we don't use it anymore
   
    ax.YGrid = 'on';        % Turn on Y-axis grid
    ax.XGrid = 'off';       % (Optional) Turn off X-axis grid
    
    ax.GridLineStyle = '--';            % Dashed grid lines
    ax.GridColor = [0.6, 0.6, 0.6];     % Light gray
    ax.GridAlpha = 0.5;                 % Semi-transparent
    
    % Force shorter xlimit
    xlim([min(time_1_i), max(time_1_i)]);
    
    % --- NOT USED ---
    % Force dynamic y-lim based on data
    % Combine all data used in the plot
    % all_vals = [conc_1_interp_i(:); conc_1_red_interp_i(:); conc_2_interp_i(:); conc_2_red_interp_i(:)];
    % all_vals = all_vals(~isinf(all_vals) & ~isnan(all_vals) & all_vals > 0);
    % all_vals = log10(all_vals);
    % Get the maximum finite value (avoid Inf/NaN)
    % max_val = max(all_vals);
    % min_val = min(all_vals);
    % Set Y limit with +20% margin
    % ylim([0.9 * min_val, 1.1 * max_val]); 
    % --- NOT USED ---
    
    %----------------------------
    savefig(fig, [char(species_red(i)), '.fig']);
end


%==========================================================================
%==========================================================================
