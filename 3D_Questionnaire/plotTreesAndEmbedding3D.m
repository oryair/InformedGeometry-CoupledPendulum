function plotTreesAndEmbedding3D(savefigs, rowtitle, coltitle, trialtitle, ...
    eigsnum_row, eigsnum_col, eigsnum_trials, ...
    row_aff, col_aff, trials_aff, ...
    row_thresh, col_thresh, trials_thresh, ...
    row_tree, col_tree, trials_tree, ...
    row_perm, col_perm, trials_perm, row_colorparam, col_colorparam, trials_colorparam)

%% Visualization
% plot the trials tree
figure;
plotTreeWithColors(trials_tree, trials_perm);
ind = strfind(trialtitle, '_');
if isempty(ind)
    nameTrial = trialtitle;
else
    nameTrial = trialtitle(1:ind(1)-1);
end
title([trialtitle 'Tree']);
if savefigs
    print(fullfile(figspath, [nameTrial 'Tree.pdf']),'-dpdf');
    saveas(gcf, fullfile(figspath, [nameTrial 'Tree.jpg']),'jpg');
end
% plot the Frequency tree
figure;
plotTreeWithColors(row_tree, (row_perm));
ind = strfind(rowtitle, '_');
if isempty(ind)
    nameRow = rowtitle;
else
    nameRow = rowtitle(1:ind(1)-1);
end
title([nameRow ' Tree']);
if savefigs
    print(fullfile(figspath, [rowtitle 'Tree.pdf']),'-dpdf');
    saveas(gcf, fullfile(figspath, [rowtitle 'Tree.jpg']),'jpg');
end
% plot the Time Frames tree
figure;
plotTreeWithColors(col_tree, col_perm);
ind = strfind(coltitle, '_');
if isempty(ind)
    nameCol = coltitle;
else
    nameCol = coltitle(1:ind(1)-1);
end
title([nameCol ' Tree']);
if savefigs
    print(fullfile(figspath, [coltitle 'Tree.pdf']),'-dpdf');
     saveas(gcf, fullfile(figspath, [coltitle 'Tree.jpg']),'jpg');
end


%% This is only for setting the thresholds
% figure;subplot(3,2,1);
% hist(row_aff(:),100);
% subplot(3,2,3);
% hist(col_aff(:),100);
% subplot(3,2,5);
% hist(trials_aff(:),100);
% subplot(3,2,2);
% [~, i] = sort(row_perm);
% imagesc(row_aff(i, i));colorbar;
% subplot(3,2,4);
% [~, i] = sort(col_perm);
% imagesc(col_aff(i, i));colorbar;
% subplot(3,2,6);
% [~, i] = sort(trials_perm);
% imagesc(trials_aff(i, i));colorbar;
% % keyboard;
% row_thresh = 0.3;
% col_thresh = 0.1;
% trials_thresh = 0.3;% 0.4
row_aff1 = threshold(row_aff, row_thresh);
col_aff1 = threshold(col_aff, col_thresh);
trials_aff1 = threshold(trials_aff, trials_thresh);
%% Get Final Embedding


[col_vecs, col_vals] = CalcEigs(col_aff1, eigsnum_col);
embedding = col_vecs*col_vals;
figure;

if isempty(col_colorparam)
    PlotEmbedding( embedding, col_perm, [nameCol  ]);
else
    subplot(2,1,1);
    PlotEmbedding( embedding, col_perm, [nameCol ]);
    subplot(2,1,2);
    plotEmbeddingWithColors(embedding, col_colorparam, [nameCol ])
end
if savefigs
    print(fullfile(figspath, ['Embedding' coltitle '.pdf']),'-dpdf');
    saveas(gcf, fullfile(figspath, ['Embedding' coltitle '.jpg']),'jpg');
end
% draw embedding for rows
[row_vecs, row_vals] = CalcEigs(row_aff1, eigsnum_row);
figure;

embedding = row_vecs*row_vals;

if isempty(row_colorparam)
    PlotEmbedding( embedding, row_perm, [nameRow ] );
else
    
    subplot(2,1,1);
    PlotEmbedding( embedding, row_perm, [nameRow ] );
    subplot(2,1,2);
    plotEmbeddingWithColors(embedding, row_colorparam, [nameRow ]);
end
if savefigs
    print(fullfile(figspath, ['Embedding' rowtitle '.pdf']),'-dpdf');
    saveas(gcf, fullfile(figspath, ['Embedding' rowtitle '.jpg']),'jpg');
end

[trials_vecs, trials_vals] = CalcEigs(trials_aff1, eigsnum_trials);
embedding = trials_vecs*trials_vals;

figure;
if isempty(trials_colorparam)
    PlotEmbedding( embedding, trials_perm, [nameTrial ] );
else
    subplot(2,1,1);
    PlotEmbedding( embedding, trials_perm, [nameTrial ] );
    subplot(2,1,2);
    plotEmbeddingWithColors(embedding, trials_colorparam, [nameTrial ])
end
if savefigs
    print(fullfile(figspath, ['Embedding' trialtitle '.pdf']),'-dpdf');
    saveas(gcf, fullfile(figspath, ['Embedding' trialtitle '.jpg']),'jpg');
end
end
