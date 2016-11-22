function [ aff_mat, slices ] = CalcInitAff( data, params, dim )
% Calculates the affinity between slices by columns of the data as thresholded
% cosine similarity between columns, or as a Gaussian kernel of width
% eps*(median distance between 5 nearest neighbors of all points).
%
% Inputs:
%     data - M-by-N-by-nT matrix
%     params - struct array with user parameters
%
% Output:
%     aff_mat - N-by-N symmetric matrix of non-negative affinities
%--------------------------------------------------------------------------
dimLen = length(size(data));
% backwards compatibility
if ~isfield(params, 'RangeMinus1to1')
    params.RangeMinus1to1 = false;
end

data   = permute(data, [setdiff(1:dimLen, dim), dim]);
slices = reshape(data, [], size(data, dimLen));
    
if dimLen >= 2
%     slices = zeros(size(data, 1)* size(data, 3), size(data, 2));
%     for r=1:size(data, 2)
%         currslice = data(:, r, :);
%         slices(:, r) = currslice(:);
%     end
    switch params.metric
        case 'cosine_similarityOnTrials'
            for T=1:size(data, 3)
                v = permute(data(:, :, T), [1 2 3]);
                inner_products = v.'*v;
                [ij_norm, ji_norm] = meshgrid(diag(inner_products));
                aff_matAll(:,:,T) = inner_products./sqrt(ij_norm.*ji_norm);
            end
            aff_mat = mean(aff_matAll, 3);
            aff_mat(aff_mat < params.thresh) = 0;
    
        case 'cosine_similarity'
        inner_products     = slices.' * slices;
        [ij_norm, ji_norm] = meshgrid( diag(inner_products) );
        aff_mat            = inner_products ./ sqrt(ij_norm .* ji_norm);
        aff_mat(aff_mat < params.thresh) = 0;
        case 'euc'
        
        euc_dist = squareform(pdist(slices.'));
        [~, nn_dist] = knnsearch(slices.', slices.', 'k', params.knn);
        sigma = params.eps * median(nn_dist(:));
        aff_mat = exp(-euc_dist.^2/(2*sigma^2));
        
        otherwise
            error('params.metric is not well define, please fix it');
            
    end
    
    if params.RangeMinus1to1  
        aff_mat = 2*(aff_mat-min(aff_mat(:)))./(max(aff_mat(:)-min(aff_mat(:))))-1;
    end
else
%     for dim = 1:size(data, 1)
%         for r=1:size(data, 3)
%             currslice = data(dim, :, r, :);
%             slices(:, r, dim) = currslice(:);
%         end
%     end
    if strcmp(params.metric,'cosine_similarity'),
        error('Not implemented yet');
        %         for d = 1
        %             inner_products(r1,r2) = currslice1.'*currslice2;
        %
        %         end
        %
        %         [ij_norm, ji_norm] = meshgrid(diag(inner_products));
        %         aff_mat = inner_products./sqrt(ij_norm.*ji_norm);
        %         aff_mat(aff_mat < params.thresh) = 0;
    else
        for r1 = 1:size(slices, 2)
            for r2 = 1:size(slices, 2)
                slice1 = permute(slices(:, r1, :), [1 3 2]);
                slice2 = permute(slices(:, r2, :), [1 3 2]);
                euc_dist(r1,r2) = sum(diag(pdist2(slice1.', slice2.')));
            end
        end
        nn_dist = sort(euc_dist.').';
        params.knn = min(params.knn, size(nn_dist, 2));
        sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
        aff_mat = exp(-euc_dist.^2/(2*sigma^2));
        
    end
end
% figure, imagesc(aff_mat), colormap gray, colorbar, axis on, hold on
% if params.on_rows,
%     title('Initial Row Affinity'), hold off
% else
%     title('Initial Col Affinity'), hold off
% end

end

