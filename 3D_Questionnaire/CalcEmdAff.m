function [ aff_mat, emd_mat ] = CalcEmdAff(data, cTrees, cParams, dim)

dim_length = ndims(data);
%-- moving the dim of interest to last:
dim_order = [setdiff(1: dim_length, dim), dim];
data      = permute(data, dim_order);
cTrees    = {cTrees{dim_order}};
cParams   = {cParams{dim_order}};

vTreesLength = zeros(1, dim_length-1);
for dim = 1 : (dim_length-1)
   vTreesLength(dim) = length(cTrees{dim});
end

emd_mat = zeros(size(data, dim_length));

%%
% Average
coefs = FindTreeAverages(data, cTrees);

P = prod(vTreesLength);
for ii = 1 : P
    
    cLevels      = cell(1, dim_length - 1);
    [cLevels{:}] = ind2sub(vTreesLength, ii);
    
    tCoefs = coefs{cLevels{:}};
    
    prod_factor = 1;
    for dim = 1 : (dim_length - 1)
        vW = cTrees{dim}{cLevels{dim}}.folder_sizes;
        vW = vW / sum(vW);
        vW = vW .^ cParams{dim}.beta;
        
        vOrder = [dim, setdiff(1: dim_length, dim)];
        tCoefs = permute(tCoefs, vOrder);
        vShape = size(tCoefs);
        mCoefs = reshape(tCoefs, vShape(1), []);
        mCoefs = bsxfun(@times, mCoefs, vW');
        tCoefs = reshape(mCoefs, vShape);
        tCoefs = ipermute(tCoefs, vOrder);
        
        prod_factor = prod_factor * 2^(cParams{dim}.alpha * (cLevels{dim} - vTreesLength(dim)));
    end
    
    w  = tCoefs;
    w2 = reshape(w, [], size(w, dim_length));
    final_W = squareform( pdist(w2', 'cityblock') );

    emd_mat = emd_mat + final_W * prod_factor;
end

eps     = cParams{dim_length}.eps * median(emd_mat(:));
aff_mat = exp(-emd_mat / eps);

end


