function [ thresholded_x ] = threshold( x, thresh )
% Thresholds a vector or matrix by setting to zero all its entries which 
% are less than threshold.
%
% Input:
%     x - either matrix or vector
%     thresh - desired threshold
% 
% Output:
%     thresholded_x - thresholded matrix or vector at the same size as x 
%--------------------------------------------------------------------------

thresholded_x = x .* (x > thresh);

end

