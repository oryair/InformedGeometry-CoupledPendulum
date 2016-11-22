function [averages, stds] = FindTreeAveragesAndStds3D( data, tree_cols, tree_rows )
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
stds{1,1} = ones(size(data));

% first average using all folders of the col tree and the bottom level of
% the rows tree - the data matrix
m = size(data, 1);
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
    mCurr = zeros(m, n, size(data,3));    
    mPrev = averages{1, col_level-1};
    stdCurr = zeros(m, n, size(data,3));  
    for ii = 1:n,
        mIdx =  tree_cols{col_level-1}.super_folders == ii;
        mCurr(:, ii, :) = mean( mPrev(:, mIdx, :), 2);
        for T=1:size(mCurr,3)
            stdCurr(:, ii, T) = std( mPrev(:, mIdx, T).').';
            mCurr1(:, ii, T) = mean( mPrev(:, mIdx, T).').';
        end
    end
    
    averages{1, col_level} = mCurr;
    stdCurr(stdCurr==0) = 1;
    stds{1, col_level} = stdCurr;
end


% first average using all folders of the row tree and the bottom level of
% the col tree - the data matrix
n = size(data, 2);
for row_level = 2:rows_n_levels,
    
    m = tree_rows{row_level}.folder_count;
    
    mCurr = zeros(m, n, size(data,3));
    stdCurr = zeros(m, n, size(data,3));
    mPrev = averages{row_level - 1, 1};
    
    for jj = 1:m,
        mIdx = tree_rows{row_level-1}.super_folders == jj;
        mCurr(jj, :, :) = mean( mPrev(mIdx, :, :), 1 );
        for T=1:size(mCurr,3)
            stdCurr(jj, :, T) = std( mPrev(mIdx, :, T) );
        end
    end
    
    averages{row_level, 1} = mCurr;
    stdCurr(stdCurr==0) = 1;
    stds{row_level, 1} = stdCurr;
end


% going over all levels of the two trees (higher than 1) and evaluating the
% averages
for row_level = 2:rows_n_levels,
    for col_level = 2:cols_n_levels,
        
        n = tree_cols{col_level}.folder_count;
        m = tree_rows{row_level}.folder_count;
        
        mCurr = zeros(m, n, size(data,3));
        stdCurr = zeros(m, n, size(data,3));
        mPrev = averages{row_level-1, col_level-1};
        
        for ii = 1:n
            mIdx = tree_cols{col_level-1}.super_folders == ii;
            vMeanPrev = mean( mPrev(:, mIdx, :), 2 );
            
            for jj = 1:m,
                mIdx2 = tree_rows{row_level-1}.super_folders == jj;                
                mCurr(jj, ii, :) = sum(vMeanPrev(mIdx2, :), 1) / nnz(mIdx2);
            end
            mat1 = mPrev(:, mIdx, :);
            for jj = 1:m,
                mIdx2 = tree_rows{row_level-1}.super_folders == jj;
                mat2 = mat1(mIdx2, :, :);
                for T=1:size(mat2,3)
                    v = mat2(:,:,T);
                    stdCurr(jj,ii,T) = std(v(:));
                end
            end
        end
        
        averages{row_level, col_level} = mCurr;
        stdCurr(stdCurr==0) = 1;
        stds{row_level, col_level} = stdCurr;
    end
    
end
end
