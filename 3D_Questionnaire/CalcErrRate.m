function [ err_rate ] = CalcErrRate( orig_data, recovered_data )
% Calculates the error rate between the original and recovered permutation 
% as the relative number of data entries that disagree
%
% Inputs: 
%     orig_data - M-by-N realization matrix  
%     recovered_data - M-by-N recovery of orig_data
%      
% Output:
%     err_rate - scalar between 0 and 1. the greater it is, the more
%         incorrect recovered permutations are
%--------------------------------------------------------------------------

err_rate = sum(sum(abs(recovered_data - orig_data)/2))/(size(orig_data,1)*size(orig_data,2));
% err_rate = pdist2(orig_data(:).', recovered_data(:).', 'hamming');
% disp(['Error rate: ',num2str(err_rate)]);

end

