function [ row_tree, col_tree, row_dual_aff, col_dual_aff ] = RunGenericQuestionnaire2D( params, data )
% Runs questionnaire algorithm using trees built by PCA clustering
% Author: Hadas Benisty
% 14.2.2016
%
% Inputs:
%     params - struct array with all user parameters
%     data - M-by-N matrix whose columns are N points in R^M (optional)
%
% Outputs:
%     row_tree, col_tree - the resulting partition trees, given as cells of
%         struct arrays
%--------------------------------------------------------------------------

if params.data.to_normalize,
    data = NormalizeData(data, params.data);
    figure, imagesc(data), colormap gray, axis on, title('Normalized Data'), colorbar
end
col_init_aff = feval(params.init_aff_col.initAffineFun, data, params.init_aff_col);

col_tree = feval(params.col_tree.buildTreeFun, col_init_aff, params.col_tree);

if params.col_tree.runOnEmbdding && params.verbose == 2
    [vecs, vals] = CalcEigs(col_init_aff, params.col_tree.eigs_num); 
    col_embedding = vecs*vals;
    figure;
    subplot(2,2,1);
    PlotEmbedding( col_embedding, 1:size(col_init_aff,1),  'Col Embedding' );
end


row_init_aff = feval(params.init_aff_row.initAffineFun, permute(data, [2 1 3]), params.init_aff_row);

row_tree = feval(params.row_tree.buildTreeFun, row_init_aff, params.row_tree);
if params.row_tree.runOnEmbdding && params.verbose == 2
    [vecs, vals] = CalcEigs(row_init_aff, params.row_tree.eigs_num); 
    row_embedding = vecs*vals;
    subplot(2,2,3);    PlotEmbedding( row_embedding, 1:size(row_init_aff,1),  'Row Embedding' );
end

if params.verbose == 2
    
    
    % plot affins
    subplot(2,2,2);    imagesc(col_init_aff);colorbar;title('Initial Col. Affin.');
    subplot(2,2,4);    imagesc(row_init_aff);colorbar;title('Initial Raw Affin.');
    % plot trees
    figure;    subplot(2,1,1);
    plotTreeWithColors(col_tree, 1:length(col_init_aff));    title('Initial Col Tree');
    subplot(2,1,2);    plotTreeWithColors(row_tree, 1:length(row_init_aff));    title('Initial Row Tree')
end


row_dual_aff = feval(params.row_tree.CalcAffFun, permute(data, [2 1 ]), col_tree, params.col_emd);
row_tree = feval(params.row_tree.buildTreeFun, row_dual_aff, params.row_tree);
figure;
for ii = 1:params.n_iters-1,
    
    col_dual_aff = feval(params.row_tree.CalcAffFun, permute(data, [1 2 ]), row_tree, params.row_emd, params.col_emd);
    col_tree = feval(params.col_tree.buildTreeFun, col_dual_aff, params.col_tree);
    if params.verbose == 2
        subplot(2,1,1);
        
        plotTreeWithColors(col_tree, 1:length(col_dual_aff));
        title(['Col Tree (Iteration ', num2str(ii),')']);
        drawnow;
    end
    row_dual_aff = feval(params.row_tree.CalcAffFun, permute(data, [2 1 ]), col_tree, params.col_emd);
    row_tree = feval(params.row_tree.buildTreeFun, row_dual_aff, params.row_tree);
    
    
    
end
if params.n_iters == 1
    col_dual_aff = feval(params.row_tree.CalcAffFun, permute(data, [1 2 ]), row_tree, params.row_emd, params.col_emd);
end

if params.verbose == 2
    % plot final affins
    figure;
    subplot(2,1,1);
    imagesc(col_dual_aff);colorbar;title('Final Col. Affin.');
    subplot(2,1,2);
    imagesc(row_dual_aff);colorbar;title('Final Raw Affin.');
    
end
end


