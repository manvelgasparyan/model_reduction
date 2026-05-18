function [] = plot_sr_graph (G, n_s, n_r)
%----------
figure
pG = plot(G, 'Layout', 'force');
pG.EdgeFontSize = 7; 
%Remove the initial labels
pG.NodeLabel = {};
%Set label positions manually
for i = 1:numel(pG.XData)
   text(pG.XData(i), pG.YData(i), G.Nodes.Name{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
%Customize the appearance of species
highlight(pG, 1:n_s, 'NodeColor', [0.0 0.6 1.0], 'Marker', 'o', 'MarkerSize', 20);
%Customize the appearance of reactions
highlight(pG, n_s+1:n_s+n_r, 'NodeColor', [1 0.8 0.0], 'Marker', 's', 'MarkerSize', 20);
%Customize the appearance of edges
set(pG, 'EdgeColor', [0, 0, 0], 'LineWidth', 1, 'ArrowSize', 5, 'EdgeAlpha', 1);
%----------
end
