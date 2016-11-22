function [ aff_mat ] = CalcEmdAff( data, tree, params, varagin) %#ok<INUSD>
% Calculates the EMD between all pairs of points (columns), and converts it 
% to an affinity.
%
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M 
%     tree - the resulting partition tree on the rows, given as cell of
%         struct arrays 
%     params - struct array with user parameters 'alpha', 'beta' and 'eps'
%     varargin - to be compatible with the generic implemetation
% 
% Output:
%     aff_mat - N-by-N symmetric matrix of non-negative affinities 
%--------------------------------------------------------------------------

N = size(data,2);
n_levels = length(tree);
emd_mat = zeros(N);

coefs = FindTreeAverages(data, tree);

for level = 1:n_levels,
    W = (tree{level}.folder_sizes / sum(tree{level}.folder_sizes)) .^ params.beta;
    W = squareform(pdist(repmat(W,N,1).*coefs{level}.', 'cityblock'));
    emd_mat = emd_mat + 2^(params.alpha*(level-n_levels)) * W;
end

eps = params.eps * median(emd_mat(:));
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

