function [ X_stoch ] = stochastic( X )
% Converts a given matrix to a right stochastic matrix.
% 
% Input:
%     X - real matrix with nonnegative entries  
% 
% Output:
%     X_stoch - right stochastic martix, with each row summing to 1  
%--------------------------------------------------------------------------

vD = sum(X, 2);
% X_stoch = D\X;
X_stoch = bsxfun(@rdivide, X, vD);

end

