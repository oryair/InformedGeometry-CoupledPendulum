function [ aff_mat, emd_mat ] = CalcEmdAff3D( data, col_tree, row_tree, params_row, params_col, params_trials, verbose )
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
T = size(data, 3);
row_n_levels = length(row_tree);
col_n_levels = length(col_tree);
emd_mat = zeros(T);


% Average
coefs  = FindTreeAverages3D(data, col_tree, row_tree );
if verbose
pb = CmdLineProgressBar('Evaluating EMD - norm 1 ');
end
for col_level = 1:col_n_levels,
    for row_level = 1:row_n_levels,
        if verbose
            pb.print(row_level + (col_level-1)*row_n_levels,row_n_levels*col_n_levels);
        end
        W_row = (row_tree{row_level}.folder_sizes / sum(row_tree{row_level}.folder_sizes)) .^ (params_row.beta);
        W_col = (col_tree{col_level}.folder_sizes / sum(col_tree{col_level}.folder_sizes)) .^ (params_col.beta);
        for T = 1:size(coefs{row_level, col_level}, 3)
            for n = 1:length(W_row)
                for m = 1:length(W_col)
                    w(n, m, T) = coefs{row_level, col_level}(n, m, T)*W_row(n)*W_col(m);
                end
            end
        end
        if any(isnan( w(:)))
            keyboard;
        end
        
%         final_W2 = zeros(size(w, 3));
%         for Ti = 1:size(w, 3)
%             for Tj = 1:size(w, 3)
%                 final_W2(Ti, Tj) = sum(sum(abs(w(:, :, Ti) - w(:, :, Tj))));
%             end
%         end
        
        w2 = reshape(w, size(w, 1) * size(w, 2), []);
        final_W = squareform( pdist(w2', 'cityblock') );
%         figure; imagesc(final_W); colorbar;
%         D = final_W - final_W2;
%         max(abs(D(:)))
        
        if any(isnan( final_W(:)))
            keyboard;
        end
        emd_mat = emd_mat + 2^(params_row.alpha*(row_level-row_n_levels)) * ...
            2^(params_col.alpha*(col_level-col_n_levels)) * final_W;
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


