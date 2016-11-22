function plotEmbeddingWithColors(embedding, colorparam_c, title_str)

color_mat     = jet(256);
colorparam_c1 = 1 + floor(255 * (colorparam_c - min(colorparam_c)) / (max(colorparam_c) - min(colorparam_c)));
if any(isnan((colorparam_c1)))
    colorparam_c1 = 1:length(colorparam_c);
end
color_mat     = color_mat(colorparam_c1, :);
switch size(embedding, 2)
    case 1
        plot(colorparam_c, embedding, 'x'), %'MarkerEdgeColor', 'k'
        title([title_str ' Embedding']);
        xlabel(title_str); ylabel('\psi_1');
%         colorbar
    case 2
        scatter(embedding(:,1), embedding(:,2), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
        title([title_str ' Embedding']);xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
        colorbar
    otherwise
        scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 60, colorparam_c, 'Fill'); colorbar;
%         scatter3(embedding(:,7), embedding(:,8), embedding(:,9), 60, colorparam_c, 'Fill'); colorbar;
        title(title_str);
        xlabel('\psi_1'), ylabel('\psi_2'); zlabel('\psi_3');
%         scatter3(embedding(:,3), embedding(:,2), embedding(:,1), 60, colorparam_c, 'Fill'); colorbar;
%         title(title_str);
%         xlabel('\psi_3'), ylabel('\psi_2'); zlabel('\psi_1');
        
end

end