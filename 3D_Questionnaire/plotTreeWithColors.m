function plotTreeWithColors(tree, labels)

nT = length(labels);

treeplot(nodes(tree), '.')
hold on;

[x_layout,y_layout] = treelayout(nodes(tree));

%-- Remove nodes and take only leafs of tree:
x_layout(1 : end - nT) = [];
y_layout(1 : end - nT) = [];

% scatter(x_layout, y_layout, 50, 1 : length(x_layout), 'fill');
scatter(x_layout, y_layout, 50, labels, 'fill');
colormap jet;
colorbar;

hold off;

end
