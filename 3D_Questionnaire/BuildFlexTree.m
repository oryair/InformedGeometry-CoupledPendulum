function [ dfs_levels,eigvals ] = BuildFlexTree( affinity, params )
% Creates a flexible tree on either the rows or the columns by calculating 
% the diffusion on the given affinity.
%
% Inputs:
%     affinity - N-by-N (symmetric) affinity matrix with non-negative entries 
%     params - struct array with user parameters 'constant' and 'min_joins_percentage'
% 
% Outputs:
%     dfs_levels - cell of struct arrays with the following entities, 
%         representing the flexible tree:
%         clustering:           vector of length N containing the folder association per level
%         super_folders:        the clustering at the next higher level of the tree
%         folder_count:         the number of folders at each level
%         folder_sizes:         vector of 'folder_count' length containing the size 
%                               of each folder, i.e. how many points are below its node
%--------------------------------------------------------------------------
if size(affinity, 1) ~= size(affinity, 2)
    error('affinity should be N by N matrix!');
end
if max(max(abs(affinity-affinity.')))> eps
    error('affinity should be symmetric');
end
if any(affinity(:) < 0)
   error('affinity should have non-negative entries'); 
end

[vecs, vals] = CalcEigs(affinity, params.eigs_num); 
eigvals      = diag(vals);
%initialize bottom-most level:
N = size(affinity,1);
dfs_levels{1}.folder_count = N;
dfs_levels{1}.folder_sizes = ones(1,N);
dfs_levels{1}.clustering   = 1:N;

dfs_time       = 1;
Q              = eye(N);
cur_level_ind  = 1;
MAX_ITERATIONS = 1e3;
isfinished     = false;

for l = 1:MAX_ITERATIONS
    if (dfs_levels{cur_level_ind}.folder_count == 1),
        dfs_levels{cur_level_ind}.super_folders = [];
        isfinished = true;
        break;
    end
    
    dfs_vecs = vecs * vals.^dfs_time;
    dfs_dists = squareform(pdist(dfs_vecs)); %Euclidean distance
    avg_dists = Q*dfs_dists*Q.'; 
    
    [clustering, redundant] = ClustByDist(avg_dists, params);
    if ~redundant,
        dfs_levels{cur_level_ind}.super_folders = clustering;
        
        cur_level_ind = cur_level_ind + 1;
        dfs_levels{cur_level_ind}.clustering = ...
            dfs_levels{cur_level_ind-1}.super_folders(dfs_levels{cur_level_ind-1}.clustering);
        dfs_levels{cur_level_ind}.folder_count = max(clustering);
        dfs_levels{cur_level_ind}.folder_sizes = CountOccurrences(clustering);
        
        Q = stochastic(ConstructClustMat(dfs_levels{cur_level_ind}.clustering));
    end
    
    dfs_time = 2*dfs_time;
end
if ~isfinished
    error('Couldn''t build tree :(');
end
% figure, treeplot(nodes(dfs_levels),'.'), title(title_str)

[~, dfs_levels] = relabel(dfs_levels);

end


function [ res ] = CountOccurrences( x )
% Counts the occurrences of each element in an indexes array.
%
% Input:
%     x - indexes vector running from 1 to D
% 
% Output:
%     res - 1-by-D vector with the number of occurrences of each index
%--------------------------------------------------------------------------

n = max(x);
res = nan(1,n);
for ii = 1:n,
    res(ii) = sum(x == ii);
end

end


function [ clust_mat ] = ConstructClustMat( clustering )
% Constructs a binary matrix of folder association.
%
% Input:
%     clustering - vector of length N, with entries in 1..F, containing the  
%         folder association, in the sense that clustering(i) is the folder
%         containing the i-th point
% 
% Output:
%     clust_mat - F-by-N binary matrix indicating the folder association, 
%         in the sense that clust_mat(j,i) equals 1 iff the i-th point 
%         belongs to the j-th folder
%--------------------------------------------------------------------------

N = length(clustering);
F = max(clustering);
clust_mat = zeros(F,N);

for ii = 1:N,
    clust_mat(clustering(ii),ii) = 1;
end

end

