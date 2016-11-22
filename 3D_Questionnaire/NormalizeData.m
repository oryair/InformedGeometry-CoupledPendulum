function [ normalized_dataT ] = NormalizeData( data, params )
% Normalizes the data over rows or columns in one of the following techniques:
%   - scaling between 0 and 1
%   - summing to unity
%   - by standard deviation
%
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M
%     params - struct array with user parameters 'normalization_dim' and 'normalization_type'
%
% Output:
%     normalized_data - M-by-N normalized matrix
%--------------------------------------------------------------------------

normalized_dataT = zeros(size(data));
for T = 1:size(data, 3)
    if ~params.over_rows,
        
        currslice = data(:,:,T).';
    else
        currslice = data(:,:,T);
    end
    switch params.normalization_type
        case 'scaling_between_0_and_1'
            normalized_data = diag(max(currslice,[],2))\currslice;
        case 'summing_to_unity'
            normalized_data = MakeStochastic(currslice);
        case 'by_std'
            normalized_data = bsxfun(@minus, currslice, mean(currslice));
            normalized_data = bsxfun(@rdivide, normalized_data, std(currslice)+eps);
        otherwise
            warning('unexpected normalization type. no normalization has been done.');
            normalized_data = currslice;
    end
    if ~params.over_rows,
        normalized_data = normalized_data.'; 
    end
    normalized_dataT(:, :, T) = normalized_data;

end

end
