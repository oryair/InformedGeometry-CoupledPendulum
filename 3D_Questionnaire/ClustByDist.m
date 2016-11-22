function [ clustering, redundant ] = ClustByDist( dist_mat, params )
% Clusters folders by diffusion distances. The resulting clusters are the
% folders at the next higher level of the tree.
%
% Inputs:
%     dist_mat - N-by-N matrix of diffusion distances between all pairs of points
%     params - struct array with user parameters 'constant' and 'min_joins_percentage'
% 
% Outputs:
%     clustering - vector of length N, with entries in 1..F, containing the  
%         folder association at this level, in the sense that clustering(i)   
%         is the folder containing the i-th point
%     redundant - boolean variable. true iff not enough joins have occurred 
%         at this level
%--------------------------------------------------------------------------

penalty = median(dist_mat(:))/params.constant; 

N = size(dist_mat,1);
dist_mat = dist_mat + realmax*eye(N);
[mins_per_row, mins_loc] = min(dist_mat,[],2);   
min_n_joins = floor(params.min_joins_percentage*N);

n_joins = 0;
clustering = 1:N;
redundant = false;

while 1,
    [min_dist, row] = min(mins_per_row);
    col = mins_loc(row);
    
    if (min_dist == realmax) && (n_joins <= min_n_joins),
        if params.verbose > 0
        disp('A redundant level - not enough joins have occurred');
        end
        redundant = true;
        break; 
    end
    if (min_dist > penalty) && (n_joins > min_n_joins),
        break; end
    
    if JoinToCluster(min_dist, clustering, row, col, penalty),
        mins_per_row(row) = realmax;
        mins_per_row(col) = realmax;
        n_joins = n_joins + 1;
        
        clustering(row) = clustering(col);
    else
        dist_mat(row,col) = realmax;
        [mins_per_row(row), mins_loc(row)] = min(dist_mat(row,:));
    end
end

[~,~,clustering] = unique(clustering);
clustering = reshape(clustering,1,[]);
        
end


function [ areJoined ] = JoinToCluster( dist, clustering, i, j, penalty )
% Attempts to join points i and j.
%
% Inputs:
%     dist - the diffusion distance between the points i and j
%     clustering - vector of length N, with entries in 1..F, containing  
%         the folder association in the sense that clustering(i) is the folder    
%         containing the i-th point 
%     i,j - indexes in the range of 1 to N
%     penalty - conatrains the size of clusters
% 
% Output:
%     areJoined - boolean variable. true iff the points have been joined
%--------------------------------------------------------------------------

Ci = sum(clustering == clustering(i));
Cj = sum(clustering == clustering(j));
max_allowable_dist = penalty*2^(1-Ci*Cj);

if (Ci+Cj == 2) || (dist <= max_allowable_dist),
    areJoined = true;
else
    areJoined = false;
end

end

