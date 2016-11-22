function newOrder = OrganizeEmbedding(data, embedding)

n_starts   = 10;
starts     = randi(min(size(data)), 1, n_starts);
l1_dist    = zeros(1, n_starts);
orders     = cell(1,  n_starts);

for ii = 1 : n_starts
    newOrder      = nn_param(embedding, starts(ii));
    
    organized_data = data(newOrder,:);
    
    [~, row_sp]    = max(sum(abs(organized_data - circshift(organized_data,[-1,0])),2));
    
    newOrder       = nn_param(embedding, newOrder(row_sp));
    organized_data = data(newOrder,:);
    
    row_sp         = sum(sum(abs(organized_data - circshift(organized_data,[-1,0])),2));
    l1_dist(ii)    = row_sp;
    orders{ii}     = newOrder;
end

[~, jj] = min(l1_dist);
newOrder = orders{jj};

end


function order = nn_param( embedding, start )
% Walks around a collapse diffusion embedding curve by taking the nearest 
% neighbor not already marked.
% 
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M
%     row_embedding, col_embedding - 
% 
% Outputs:
%     row_order, col_order -  
%--------------------------------------------------------------------------

N                 = size(embedding,1);
neighbors_indexes = knnsearch(embedding, embedding, 'k', N);
order             = start; 

for ii = 1:N,
    temp = neighbors_indexes(order(end),:);
    nn = temp(~ismember(temp, order));
    if ~isempty(nn),
        order(end+1) = nn(1);
    end
end

end

