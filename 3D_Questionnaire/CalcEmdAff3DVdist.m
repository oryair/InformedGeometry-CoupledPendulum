function [ aff_mat ] = CalcEmdAff3DVdist( data, col_tree, row_tree, params_row, params_col, params_trials, on_rows, iter )
% Calculates the EMD between all pairs of points (defined by the third dim.), and converts it
% to an affinity.
%
% Inputs:
%     data - M-by-N-by-T matrix
%     col_tree, row_tree - the resulting partition tree on the cols/rows, given as cell of
%         struct arrays
%     params - struct array with user parameters 'alpha', 'beta' and 'eps'
%     on_rows - boolean variable. if true, the affinity is on the rows (optional)
%     iter - iteration number (optional)
%
% Output:
%     aff_mat - T-by-T symmetric matrix of non-negative affinities
%--------------------------------------------------------------------------
dimLen = length(size(data));
T = size(data,dimLen);
row_n_levels = length(row_tree);
col_n_levels = length(col_tree);
emd_mat = zeros(T);


% Average
coefs  = FindTreeAverages3D(data, col_tree, row_tree );
M=[];
for col_level = 2:col_n_levels,
    for row_level = 2:row_n_levels,
        
        
        W_row = repmat((row_tree{row_level}.folder_sizes / sum(row_tree{row_level}.folder_sizes)) .^ (params_row.beta), size(data,1), 1);
        W_col = repmat((col_tree{col_level}.folder_sizes / sum(col_tree{col_level}.folder_sizes)) .^ (params_col.beta), size(data,1), 1);
        for T = 1:size(coefs{row_level, col_level}, 4)
            for n = 1:size(W_row, 2)
                for m = 1:size(W_col, 2)
                    w(:, n, m, T) = coefs{row_level, col_level}(:, n, m, T).*W_row(:, n).*W_col(:, m);
                end
            end
        end
        if any(isnan( w(:)))
            keyboard;
        end
        
        final_W = zeros(size(w, 1), size(w, 4),size(w, 4));
        for Ti = 1:size(w, 4)
            for Tj = 1:size(w, 4)
                final_W(:, Ti, Tj) = sum(sum(abs(w(:, :, :, Ti) - w(:, :, :, Tj)),2),3);
            end
        end
        
       
        if any(isnan( final_W(:)))
            keyboard;
        end
        % for now, summing over all state-space dims
        M = cat(3, M, 2^(params_row.alpha*(row_level-row_n_levels)) * ...
            2^(params_col.alpha*(col_level-col_n_levels)) * shiftdim(sum(final_W, 1)));
        emd_mat = emd_mat + 2^(params_row.alpha*(row_level-row_n_levels)) * ...
            2^(params_col.alpha*(col_level-col_n_levels)) * shiftdim(sum(final_W, 1));
        if any(isnan( emd_mat(:)))
            keyboard;
        end
    end
end
eps = params_trials.eps * median(emd_mat(:));

aff_mat = exp(-emd_mat/eps);

% if nargin > 3,
%     figure, imagesc(aff_mat), colormap gray, colorbar, axis on, hold on
%     if on_rows,
%         title(['EMD Row Affinity (Iteration ',num2str(iter),')']), hold off
%     else
%         title(['EMD Col Affinity (Iteration ',num2str(iter),')']), hold off
%     end
% end

end

function [ averages ] = FindTreeAverages( data, tree )
% Averages each data point across different scales (tree levels).
%
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M
%     tree - the resulting partition tree on the rows, given as cell of
%         struct arrays
%
% Output:
%     averages - struct array with the mean values of each point on each folder
%--------------------------------------------------------------------------

data(isnan(data)) = 0;
n_levels = length(tree);
averages = cell(1,n_levels);
averages{1} = data;

for level = 2:n_levels,
    n = tree{level}.folder_count;
    averages{level} = zeros(n, size(data,2));
    
    for ii = 1:n,
        averages{level}(ii,:) = mean(averages{level-1}((tree{level-1}.super_folders == ii),:));
    end
end

end

function averages = FindTreeAverages3D( data, tree_cols, tree_rows )
% Averages each data point across different scales (tree levels).
%
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M
%     tree - the resulting partition tree on the rows, given as cell of
%         struct arrays
%
% Output:
%     averages - struct array with the mean values of each point on each folder
%--------------------------------------------------------------------------

data(isnan(data)) = 0;
cols_n_levels = length(tree_cols);
rows_n_levels = length(tree_rows);

averages = cell(rows_n_levels,cols_n_levels);
averages{1,1} = data;

% first average using all folders of the col tree and the bottom level of
% the rows tree - the data matrix
m = size(data, 2);
% for col_level = 2:cols_n_levels,
%     
%     n = tree_cols{col_level}.folder_count;
%     averages{1, col_level} = zeros(m, n, size(data,3));
%     
%     for ii = 1:n,
%         for jj = 1:m,
%             averages{1, col_level}(jj, ii, :) = mean( mean( averages{1, col_level-1}(jj, tree_cols{col_level-1}.super_folders == ii, :), 1 ), 2 );
%         end
%     end
% end


for col_level = 2:cols_n_levels,
    
    n     = tree_cols{col_level}.folder_count;
    mCurr = zeros(size(data, 1), m, n, size(data,4));    
    mPrev = averages{1, col_level-1};
    
    for ii = 1:n,
        mIdx =  tree_cols{col_level-1}.super_folders == ii;
        
        mCurr(:, :, ii, :) = mean( mPrev(:, :, mIdx, :), 3);
    end
    
    averages{1, col_level} = mCurr;
end


% first average using all folders of the row tree and the bottom level of
% the col tree - the data matrix
n = size(data, 3);
for row_level = 2:rows_n_levels,
    
    m = tree_rows{row_level}.folder_count;
    
    mCurr = zeros(size(data,1), m, n, size(data,4));
    mPrev = averages{row_level - 1, 1};
    
    for jj = 1:m,
        mIdx = tree_rows{row_level-1}.super_folders == jj;
        mCurr(:, jj, :, :) = mean( mPrev(:, mIdx, :, :), 2 );
    end
    
    averages{row_level, 1} = mCurr;
end


% going over all levels of the two trees (higher than 1) and evaluating the
% averages
for row_level = 2:rows_n_levels,
    for col_level = 2:cols_n_levels,
        
        n = tree_cols{col_level}.folder_count;
        m = tree_rows{row_level}.folder_count;
        
        mCurr = zeros(size(data,1), m, n, size(data,4));
        mPrev = averages{row_level-1, col_level-1};
        
        for ii = 1:n
            mIdx = tree_cols{col_level-1}.super_folders == ii;
            vMeanPrev = mean( mPrev(:, :, mIdx, :), 3 );
            for jj = 1:m,
                mIdx2 = tree_rows{row_level-1}.super_folders == jj;
                mCurr(:, jj, ii, :) = sum(vMeanPrev(:, mIdx2, :), 2) / nnz(mIdx2);
            end
        end
        
        averages{row_level, col_level} = mCurr;
    end
    
end
end

