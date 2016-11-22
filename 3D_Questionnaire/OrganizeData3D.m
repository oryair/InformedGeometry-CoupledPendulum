function [ err_rate, organized_data, row_order, col_order, trial_order ] = OrganizeData3D( orig_data, data, row_aff, col_aff, trial_aff, row_perm, col_perm, trial_perm, eigsnum_col, eigsnum_row, eigsnum_trial )
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
figure;
subplot(3,1,1);
[row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
PlotEmbedding(row_vecs*row_vals, row_perm, 'Row Embedding');
subplot(3,1,2);
[col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
PlotEmbedding(col_vecs*col_vals, col_perm, 'Col Embedding');
subplot(3,1,3);
[trial_vecs, trial_vals] = CalcEigs(trial_aff, eigsnum_trial);
PlotEmbedding(trial_vecs*trial_vals, trial_perm, 'Trial Embedding');


% [~, row_order] = sort(row_vecs(:,1)); row_order = row_order.';
% [~, col_order] = sort(col_vecs(:,1)); col_order = col_order.';
% [ row_order ] = OrganizeDiffusion3DbyOneDim( data, row_vecs*row_vals );
[row_order, col_order, trial_order] = OrganizeDiffusion3D(data, row_vecs*row_vals, col_vecs*col_vals, trial_vecs*trial_vals); %rather computationally expensive 
% [row_vecs, row_vals] = CalcEigs(row_aff, 3);
% PlotEmbeddingWithColors(row_vecs*row_vals, row_order, 'Row Embedding');

[row_order, col_order, trial_order, err_rate] = CorrectOrientation3D(row_order, col_order, trial_order, orig_data, data);
organized_data = data(row_order,:, :); 
organized_data = organized_data(:,col_order, :); 
organized_data = organized_data(:,:, trial_order); 

figure, subplot(2,1,1),imagesc(orig_data(:,:,1)), colormap jet, axis on, 
title(['Orig. Data']), colorbar
subplot(2,1,2), imagesc(organized_data(:,:,1)),colormap jet, axis on, 
title(['Recovered Permutation']), colorbar;
end


function [ correct_row_order, correct_col_order, correct_trial_order, correct_err_rate ] = ...
    CorrectOrientation3D( row_order, col_order, trial_order, orig_data, data )
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
l=1;
for flip1 = [0 1]
    for flip2 = [0 1]
        for flip3 = [0 1]
            if flip1
                row_orderused = fliplr(row_order);
            else
                row_orderused = (row_order);
            end
            if flip2
                col_orderused = fliplr(col_order);
            else
                col_orderused = (col_order);
            end
            if flip3
                trial_orderused = fliplr(trial_order);
            else
                trial_orderused = (trial_order);
            end
            
            organized_data = data(row_orderused,:, :);
            organized_data = organized_data(:,col_orderused, :);
            organized_data = organized_data(:,:, trial_orderused);
            err_rate(l) = CalcErrRate(orig_data, organized_data);
            flipused(1, l)= flip1;
            flipused(2, l)= flip2;
            flipused(3, l)= flip3;
            l = l + 1;
            
        end
    end
end

[correct_err_rate, ind] = min(err_rate);
if flipused(1, ind)
    correct_row_order = fliplr(row_order);
else
    correct_row_order = row_order;
end
if flipused(2, ind)
    correct_col_order = fliplr(col_order);
else
    correct_col_order = col_order;
end
if flipused(2, ind)
    correct_trial_order = fliplr(trial_order);
else
    correct_trial_order = trial_order;
end


end
