function [ vecs, vals ] = CalcEigs( A, n_eigs )
% Displays the 3-D diffusion embedding, i.e. the 3 eigenvectors corresponding 
% to the 3 largest eigenvalues.
%
% Inputs:
%     A - N-by-N (symmetric) affinity matrix with non-negative entries 
%     n_eigs - number of desired eigenvalues/eigenvectors
% 
% Outputs:
%     vals - diagonal matrix of nontrivial eigenvalues (all smaller than 1),
%         sorted in a descend order 
%     vecs - the corresponding eigenvectors  
%--------------------------------------------------------------------------

% [stoch_A, vD] = stochastic(A);
n_eigs = min(n_eigs, size(A,1));
% [vecs, vals] = eigs(stoch_A, eigs_n); 
% [vecs, vals, ~] = svds(stoch_A, n_eigs);

vD = sum(A,2);
one_over_D_sqrt = spdiags(sqrt(1./vD),0,size(A,1),size(A,2));
% using symmetric matrix for calculation
A_sym = one_over_D_sqrt * A * one_over_D_sqrt;
A_sym = 0.5*(A_sym + A_sym'); % to force A_sym to be numerically symmetric

[vecs, vals] = eig(A_sym);
[vals, I] = sort(diag(vals), 'descend'); 
I         = I(1 : n_eigs);
vals      = vals(1 : n_eigs);
vecs      = vecs(:,I);
vecs      = one_over_D_sqrt * vecs; %convert to eigenvectors of stoch_A
% figure, plot(vals, 'r+'), title('Affinity eigenvalues')

%assuming that the induced graph is connected, the first eigenvalue equals
%1 (others decay from there) and the corresponding eigenvector is constant.
%we normalize all the eigenvectors such that the latter is 1:
v = vecs(:,1);
v(abs(v) <= 1e-7) = 1e-6;
vecs = bsxfun(@rdivide, vecs, v);

sign_vec = sign(vecs(1,:));
sign_vec(sign_vec == 0) = 1;
vecs = bsxfun(@rdivide, vecs, sign_vec);

%now we omit the trivial terms:
vals = diag(vals(2:end));
vecs = vecs(:,2:end);

end

