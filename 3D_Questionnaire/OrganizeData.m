function [ err_rate, organized_data_row_col, row_order, col_order ] = OrganizeData( orig_data, data, row_aff, col_aff, row_perm, col_perm, eigsnum_col, eigsnum_row  )
% Attemps to recover the row and column permutations by organizing the data 
% according to diffusion maps embedding (its dimensionality is set to 1, 
% thus it associates each image with a scalar).
%
% Inputs:
%     orig_data - M-by-N realization matrix  
%     data - M-by-N permutation of data
%     row_aff, col_aff - M-by-M and N-by-N dual affinity matrices
%     row_perm, col_perm - permutations on the rows and columns of the original data
% 
% Outputs:
%     err_rate - scalar between 0 and 1. the greater it is, the more
%         incorrect recovered permutations are
%--------------------------------------------------------------------------

[row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
figure;
subplot(2,1,1);
PlotEmbedding(row_vecs*row_vals, row_perm, 'Row Embedding');

[col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
subplot(2,1,2);
PlotEmbedding(col_vecs*col_vals, col_perm, 'Col Embedding');

% [~, row_order] = sort(row_vecs(:,1)); row_order = row_order.';
% [~, col_order] = sort(col_vecs(:,1)); col_order = col_order.';
[row_order, col_order] = OrganizeDiffusion(data, row_vecs*row_vals, col_vecs*col_vals); %rather computationally expensive 

[row_order, col_order, err_rate] = CorrectOrientation(row_order, col_order, orig_data, data);
organized_data_row = data(row_order,:); 
organized_data_row_col = organized_data_row(:,col_order); 
organized_data_col = data(:,col_order); 
figure, subplot(2,2,1),imagesc(orig_data), colormap jet, axis on, 
title(['Orig. Data']), colorbar
subplot(2,2,2), imagesc(organized_data_row),colormap jet, axis on, 
title(['Recovered Permutation By Row']), colorbar;
subplot(2,2,3), imagesc(organized_data_col),colormap jet, axis on, 
title(['Recovered Permutation By Col']), colorbar;

subplot(2,2,4), imagesc(organized_data_row_col),colormap jet, axis on, 
title(['Recovered Permutation By Row+Col']), colorbar;

end


function [ correct_row_order, correct_col_order, correct_err_rate ] = ...
    CorrectOrientation( row_order, col_order, orig_data, data )
% Reverse the row and column orders if needed, otherwise the recovered 
% permutation would be a mirror image of the original.
%
% Inputs:
%     row_order, col_order - resulted orders, might be flipped
%     orig_data - M-by-N realization matrix  
%     data - M-by-N permutation of data
% 
% Outputs:
%     correct_row_order, correct_col_order - properly oriented orders
%     err_rate - scalar between 0 and 1. the greater it is, the more
%         incorrect recovered permutations are
%--------------------------------------------------------------------------
organized_data = data(row_order,:); 
organized_data = organized_data(:,col_order); 
err_rate(1) = CalcErrRate(orig_data, organized_data);
  
organized_data = data(fliplr(row_order),:); 
organized_data = organized_data(:,col_order); 
err_rate(2) = CalcErrRate(orig_data, organized_data);
 
organized_data = data(row_order,:); 
organized_data = organized_data(:,fliplr(col_order)); 
err_rate(3) = CalcErrRate(orig_data, organized_data);
 
organized_data = data(fliplr(row_order),:); 
organized_data = organized_data(:,fliplr(col_order)); 
err_rate(4) = CalcErrRate(orig_data, organized_data);
     
[correct_err_rate, ind] = min(err_rate);

switch ind,
     case 1
        correct_row_order = row_order;
        correct_col_order = col_order;
    case 2
        correct_row_order = fliplr(row_order);
        correct_col_order = col_order;
    case 3
        correct_row_order = row_order;
        correct_col_order = fliplr(col_order);
    case 4
        correct_row_order = fliplr(row_order);
        correct_col_order = fliplr(col_order);
end

end
