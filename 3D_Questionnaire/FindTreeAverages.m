function cAverages = FindTreeAverages(data, cTrees)

dim_length   = ndims(data);
vTreesLength = zeros(1, dim_length-1);
N            = 1;
for dim = 1 : (dim_length-1)
   vTreesLength(dim) = length(cTrees{dim});
   N                 = N * size(data, dim);
end

cAverages    = cell(vTreesLength);
cAverages{1} = data;
cSizes       = cell(vTreesLength);
cSizes{1}    = ones(size(cAverages{1}));

%%
P       = prod(vTreesLength);
CellAdd = @(c, n) num2cell(max(cell2mat(c) + n, 1));
for ii = 2 : P
    cLevels      = cell(1, dim_length - 1);
    [cLevels{:}] = ind2sub(vTreesLength, ii);

    cIdx_prev = CellAdd(cLevels, -1);
    tData     = cAverages{cIdx_prev{:}};
    tSize     = cSizes{cIdx_prev{:}};
    tData     = tData .* tSize;
    for kk = 1 : (dim_length - 1)
        if cLevels{kk} == 1
            continue;
        end
        [tData, tSize] = TreeTensorMean(tData, tSize, cTrees{kk}, cLevels{kk}, kk);
    end
    cSizes{cLevels{:}}    = tSize;
    cAverages{cLevels{:}} = tData ./ tSize;
    
end

end
