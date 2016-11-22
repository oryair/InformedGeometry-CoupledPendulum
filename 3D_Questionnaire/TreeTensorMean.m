function [tMeandData, tMeanSize] = TreeTensorMean(tData, tSize, tree, level, dim)

    dim_length = ndims(tData);
    dim_order  = [dim, setdiff(1: dim_length, dim)];
    vDataSize  = size(tData);
    mData      = reshape(permute(tData, dim_order), vDataSize(dim), []);
    mSize      = reshape(permute(tSize, dim_order), vDataSize(dim), []);
    
    folder_count   = tree{level}.folder_count;
    vDataSize(dim) = folder_count;
    mDataCurr      = zeros(folder_count, size(mData, 2));
    mSizeCurr      = zeros(folder_count, size(mSize, 2));
    
    for folder = 1 : folder_count
        mIdx             = tree{level-1}.super_folders == folder;
        mDataCurr(folder, :) = sum(mData(mIdx, :), 1);
        mSizeCurr(folder, :) = sum(mSize(mIdx, :), 1);
    end
    
    vShape     = vDataSize(dim_order);
    vShape(1)  = folder_count;
    
    mTemp     = reshape(mDataCurr, vShape);
    tMeandData = ipermute(mTemp, dim_order);
    
    mTemp     = reshape(mSizeCurr, vShape);
    tMeanSize = ipermute(mTemp, dim_order);
    
end