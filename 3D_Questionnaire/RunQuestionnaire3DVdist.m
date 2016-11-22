function [ row_tree, col_tree, trials_tree, row_dual_aff_final, col_dual_aff_final, trials_dual_aff_final, col_init_aff, row_init_aff, trials_init_aff ] = RunQuestionnaire3DVdist( params, data )
% Runs questionnaire algorithm as suggested in [J. I. Ankenman, PhD thesis,
% Yale University, 2014].
%
% Inputs:
%     params - struct array with all user parameters
%     data - M-by-N matrix whose columns are N points in R^M (optional)
%
% Outputs:
%     row_tree, col_tree - the resulting partition trees, given as cells of
%         struct arrays
%     row_dual_aff_final, col_dual_aff_final, trials_dual_aff_final  -
%     the final affinity matrices evaluated using the EMD (using the 3D-trees-metric)
%--------------------------------------------------------------------------

% normalize 
if params.data.to_normalize,
    data = NormalizeData(data, params.data);
end


% get initial affin. matrices
[col_init_aff] = CalcInitAff3D(data, params.init_aff_col);
[col_tree, col_eigs(:,1)] = BuildFlexTree(col_init_aff, params.col_tree, 'Initial Col Tree');

[row_init_aff] = CalcInitAff3D(permute(data, [1 3  2 4]), params.init_aff_row);
[row_tree, row_eigs(:,1)] = BuildFlexTree(row_init_aff, params.row_tree, ['Initial Row Tree']);

[trials_init_aff] = CalcInitAff3D(permute(data, [1 2 4 3]), params.init_aff_trials);
[trials_tree, trials_eigs(:,1)] = BuildFlexTree(trials_init_aff, params.trials_tree, ['Initial Trails Tree']);
if params.verbose == 2
    figure;
    subplot(3,2,1);
    [vecs, vals] = CalcEigs(row_init_aff, 4);
    embedding = vecs*vals;
    PlotEmbedding( embedding, 1:size(row_init_aff,1), [ 'Row Embedding'] );
    subplot(3,2,3);
    [vecs, vals] = CalcEigs(col_init_aff, 4);
    embedding = vecs*vals;
    PlotEmbedding( embedding, 1:size(col_init_aff,1), [ 'Col Embedding'] );
    subplot(3,2,5);
    [vecs, vals] = CalcEigs(trials_init_aff, 4);
    embedding = vecs*vals;
    PlotEmbedding( embedding, 1:size(trials_init_aff,1), [ 'Trial Embedding'] );
    
    
    % plot affins
    subplot(3,2,2);
    imagesc(col_init_aff);colorbar;title('Initial Col. Affin.');
    subplot(3,2,4);
    imagesc(row_init_aff);colorbar;title('Initial Raw Affin.');
    subplot(3,2,6);
    imagesc(trials_init_aff);colorbar;title('Initial Trials Affin.');
    
    % plot trees
    figure;
    subplot(3,1,1);
    plotTreeWithColors(col_tree, 1:length(col_init_aff));
    title('Initial Col Tree');
    
    subplot(3,1,2);
    plotTreeWithColors(row_tree, 1:length(row_init_aff));
    title('Initial Row Tree')
    
    subplot(3,1,3);
    plotTreeWithColors(trials_tree, 1:length(trials_init_aff));
    title('Initial Trials Tree');
    
    
    
end

% main loop 
if params.init_aff_col.on_rows,
    error('Not implemented yet');
    for ii = 1:params.n_iters,
        trials_dual_aff = CalcEmdAff3D(permute(data, [2 1 3]), row_tree, col_tree, params.col_emd, params.row_emd, params.trials_emd, ~params.init_aff_col.on_rows, ii);
        [trials_tree, trials_eigs(:,ii+1)] = BuildFlexTree(trials_dual_aff, params.trials_tree, ['Trials Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,3);
            plotTreeWithColors(trials_tree, 1:length(trials_dual_aff));
            title(['Trials Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        
        col_dual_aff = CalcEmdAff3D(permute(data, [1 3 2]), trials_tree, row_tree, params.row_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
        [col_tree, col_eigs(:,ii+1)] = BuildFlexTree(col_dual_aff, params.col_tree, ['Col Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,1);
            
            plotTreeWithColors(col_tree, 1:length(col_dual_aff));
            title(['Col Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        
        row_dual_aff = CalcEmdAff3D(permute(data, [2 3 1]), trials_tree, col_tree, params.col_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
        [row_tree, row_eigs(:,ii+1)] = BuildFlexTree(row_dual_aff, params.row_tree, ['Row Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,2);            
            plotTreeWithColors(row_tree, 1:length(row_dual_aff));
            title(['Row Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        
        if mean(abs(diff(trials_eigs(:,end-1:end)')'./trials_eigs(:, end))) <= params.trials_tree.eigs_diff_th || ...
                mean(abs(diff(row_eigs(:,end-1:end)')'./row_eigs(:, end))) <= params.row_tree.eigs_diff_th || ...
                mean(abs(diff(col_eigs(:,end-1:end)')'./col_eigs(:, end))) <= params.col_tree.eigs_diff_th
            break;
        end
        
    end
    trials_dual_aff_final = CalcEmdAff3D(permute(data, [2 1 3]), row_tree, col_tree, params.col_emd, params.row_emd, params.trials_emd, ~params.init_aff_col.on_rows, ii);
    
else
    for ii = 1:params.n_iters,
        
        trials_dual_aff = CalcEmdAff3DVdist(data, col_tree, row_tree, params.row_emd, params.col_emd, params.trials_emd, ~params.init_aff_col.on_rows, ii);
        [trials_tree, trials_eigs(:,ii+1)] = BuildFlexTree(trials_dual_aff, params.trials_tree, ['Trials Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,3);  
            plotTreeWithColors(trials_tree, 1:length(trials_dual_aff));
            title(['Trials Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        
        row_dual_aff = CalcEmdAff3DVdist(permute(data, [1 3 4 2]), trials_tree, col_tree, params.col_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
        [row_tree, row_eigs(:,ii+1)] = BuildFlexTree(row_dual_aff, params.row_tree, ['Row Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,2);            
            plotTreeWithColors(row_tree, 1:length(row_dual_aff));
            title(['Row Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        col_dual_aff = CalcEmdAff3DVdist(permute(data, [1 2 4 3]), trials_tree, row_tree, params.row_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
        [col_tree, col_eigs(:,ii+1)] = BuildFlexTree(col_dual_aff, params.col_tree, ['Col Tree (Iteration ', num2str(ii),')']);
        if params.verbose == 2
            subplot(3,1,1);
            
            plotTreeWithColors(col_tree, 1:length(col_dual_aff));
            title(['Col Tree (Iteration ', num2str(ii),')']);
            drawnow;
        end
        if mean(abs(diff(trials_eigs(:,end-1:end)')'./trials_eigs(:, end))) <= params.trials_tree.eigs_diff_th || ...
                mean(abs(diff(row_eigs(:,end-1:end)')'./row_eigs(:, end))) <= params.row_tree.eigs_diff_th || ...
                mean(abs(diff(col_eigs(:,end-1:end)')'./col_eigs(:, end))) <= params.col_tree.eigs_diff_th
            break;
        end
    end
    trials_dual_aff_final = CalcEmdAff3DVdist(data, col_tree, row_tree, params.row_emd, params.col_emd, params.trials_emd, ~params.init_aff_col.on_rows, ii);
end

row_dual_aff_final = CalcEmdAff3DVdist(permute(data, [1 3 4 2]), trials_tree, col_tree, params.col_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
col_dual_aff_final = CalcEmdAff3DVdist(permute(data, [1 2 4 3]), trials_tree, row_tree, params.row_emd, params.trials_emd, params.col_emd, ~params.init_aff_col.on_rows, ii);
if params.verbose == 2
    % plot final affins
    figure;
    subplot(3,1,1);
    imagesc(col_dual_aff_final);colorbar;title('Final Col. Affin.');
    subplot(3,1,2);
    imagesc(row_dual_aff_final);colorbar;title('Final Raw Affin.');
    subplot(3,1,3);
    imagesc(trials_dual_aff_final);colorbar;title('Final Trials Affin.');
end