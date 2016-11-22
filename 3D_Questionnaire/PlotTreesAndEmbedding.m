function PlotTreesAndEmbedding(cAff, cTrees, cColor, cTitles)



%% Visualization
for ii = 1 : length(cAff)

    if nargin < 4
        cTitles{ii} = ['Dim ', num2str(ii)];
    end
    
    figure;
    plotTreeWithColors(cTrees{ii}, cColor{ii});
    title(cTitles{ii});
    
%     [vecs, vals] = CalcEigs(cAff{ii}, 4);
    [vecs, vals] = CalcEigs(cAff{ii}, inf);
    embedding    = vecs * vals;
    
    figure;
    if nargin < 3
        PlotEmbedding( embedding, 1 : size(vecs, 1), cTitles{ii});
    else
        plotEmbeddingWithColors(embedding, cColor{ii}, cTitles{ii})
    end

end

end
