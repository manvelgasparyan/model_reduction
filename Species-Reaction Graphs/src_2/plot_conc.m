function [] = plot_conc (time, conc, species)
%---------------
n_s = size(conc{1},1);
%---------------
for j = 1:n_s
    figure
    for i = 1:numel(conc)
        %------------------
        conc_i = conc{i};
        %------------------
        x_j_i = conc_i(j,:);
        %------------------
        plot(time, x_j_i, 'r-', 'LineWidth', 1.5); % Blue line, thickness = 2
        hold on
        %-------------------
    end
    %---
    hold off;
    title([char(species(j)) ': Original'], 'FontSize', 20, 'FontWeight', 'bold');
    %---
    xlabel('Time (seconds)', 'FontSize', 18);
    ylabel('Concentration (mol)', 'FontSize', 18);
    %---
    %%legend({'Original'}, 'FontSize', 18, 'Location', 'best');
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
%----------------
end
