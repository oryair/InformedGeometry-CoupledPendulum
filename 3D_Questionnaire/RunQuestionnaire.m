function [ row_tree, col_tree, row_dual_aff, col_dual_aff ] = RunQuestionnaire( params, data )
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
%--------------------------------------------------------------------------

if nargin < 2,
    load(params.data.file);
    figure, imagesc(data), colormap gray, axis on, title('Original Data'), colorbar
end

if params.data.to_normalize,
    data = NormalizeData(data, params.data);
    figure, imagesc(data), colormap gray, axis on, title('Normalized Data'), colorbar
end

if params.init_aff_row.on_rows,
    row_init_aff = CalcInitAff(data.', params.init_aff_row);
    row_tree = BuildFlexTree(row_init_aff, params.row_tree, 'Initial Row Tree');
    
    for ii = 1:params.n_iters,
        col_dual_aff = CalcEmdAff(data, row_tree, params.col_emd, ~params.init_aff_col.on_rows, ii);
        col_tree = BuildFlexTree(col_dual_aff, params.col_tree, ['Col Tree (Iteration ',num2str(ii),')']);
        
        row_dual_aff = CalcEmdAff(data.', col_tree, params.row_emd, params.init_aff_col.on_rows, ii);
        row_tree = BuildFlexTree(row_dual_aff, params.row_tree, ['Row Tree (Iteration ', num2str(ii),')']);
    end
else
    col_init_aff = CalcInitAff2D(data, params.init_aff_col);
    col_tree = BuildFlexTree(col_init_aff, params.col_tree);
    figure;
    
    for ii = 1:params.n_iters,
        row_dual_aff = CalcEmdAff2D(data.', col_tree, params.row_emd);
        row_tree = BuildFlexTree(row_dual_aff, params.row_tree);
        
        col_dual_aff = CalcEmdAff2D(data, row_tree, params.col_emd);
        col_tree = BuildFlexTree(col_dual_aff, params.col_tree);
        subplot(2,1,1);
        
        plotTreeWithColors(row_tree, 1:length(row_dual_aff));
        title('2D - Row Tree');
        
        subplot(2,1,2);
        plotTreeWithColors(col_tree, 1:length(col_dual_aff));
        title('2D - Col Tree');
        drawnow;
        
    end
end

end

