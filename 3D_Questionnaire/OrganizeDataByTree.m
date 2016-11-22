function [organized_data_row_col, row_order, col_order ] = OrganizeDataByTree(data, row_tree, col_tree)
% Organizes data by tree structure
%
% Inputs:
%     data - M-by-N permutation of data
%     row_tree, col_tree - trees on rows and cols
% 
% Outputs:
%     organized_data_row_col - data with rows and cols permuted by tree
%     structure
%--------------------------------------------------------------------------
[~, row_order]  = sort(row_tree{2}.clustering);
[~, col_order] = sort(col_tree{2}.clustering);

organized_data_row = data(row_order,:); 
organized_data_row_col = organized_data_row(:,col_order); 
organized_data_col = data(:,col_order); 
figure, subplot(2,2,1),imagesc(data), colormap jet, axis on, 
title(['Orig. Data']), colorbar
subplot(2,2,2), imagesc(organized_data_row),colormap jet, axis on, 
title(['Permutation By Row Tree']), colorbar;
subplot(2,2,3), imagesc(organized_data_col),colormap jet, axis on, 
title(['Permutation By Col Tree']), colorbar;

subplot(2,2,4), imagesc(organized_data_row_col),colormap jet, axis on, 
title(['Permutation By Row+Col Trees']), colorbar;

end

