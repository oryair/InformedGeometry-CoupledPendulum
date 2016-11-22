function [ row_order, col_order, trial_order ] = OrganizeDiffusion3D( data, row_embedding, col_embedding, trial_embedding )
% Recovers the original permutation of the shuffled data based on the diffusion embeddings of rows and columns.
% 
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M
%     row_embedding, col_embedding - 
% 
% Outputs:
%     row_order, col_order - 
%--------------------------------------------------------------------------

n_starts = 10;
starts = randi(min(size(data)),1,n_starts);
l1_dist = zeros(1,n_starts);
row_orders = cell(1,n_starts);
col_orders = cell(1,n_starts);
trial_orders = cell(1,n_starts);

for ii = 1:n_starts,
    row_order = nn_param(row_embedding, starts(ii));
    col_order = nn_param(col_embedding, starts(ii));
    trial_order = nn_param(trial_embedding, starts(ii));
    organized_data = data(row_order,:, :);
    organized_data = organized_data(:,col_order, :);
    organized_data = organized_data(:,:, trial_order);
    
    [~, row_sp] = max(sum(sum(abs(organized_data - circshift(organized_data,[-1,0,0])),2),3));
    [~, col_sp] = max(sum(sum(abs(organized_data - circshift(organized_data,[0,-1,0])),1),3));
    [~, trial_sp] = max(sum(sum(abs(organized_data - circshift(organized_data,[0,-1,0])),1),2));
    
    row_order = nn_param(row_embedding, row_order(row_sp));
    col_order = nn_param(col_embedding, col_order(col_sp));
    trial_order = nn_param(trial_embedding, trial_order(trial_sp));
    organized_data = data(row_order,:, :);
    organized_data = organized_data(:,col_order, :);
    organized_data = organized_data(:,:, trial_order);
    
    row_sp = sum(sum(sum(abs(organized_data - circshift(organized_data,[-1,0,0])))));
    col_sp = sum(sum(sum(abs(organized_data - circshift(organized_data,[0,-1,0])))));
    trial_sp = sum(sum(sum(abs(organized_data - circshift(organized_data,[0,0,-1])))));
    
    l1_dist(ii) = row_sp + col_sp + trial_sp;
    row_orders{ii} = row_order;
    col_orders{ii} = col_order;
    trial_orders{ii} = trial_order;
end

[~, jj] = min(l1_dist);
row_order = row_orders{jj};
col_order = col_orders{jj};
trial_order = trial_orders{jj};
end


function [ order ] = nn_param( embedding, start )
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

N = size(embedding,1);
neighbors_indexes = knnsearch(embedding, embedding, 'k', N);
order = start; 

for ii = 1:N,
    temp = neighbors_indexes(order(end),:);
    nn = temp(~ismember(temp, order));
    if ~isempty(nn),
        order(end+1) = nn(1);
    end
end

end

