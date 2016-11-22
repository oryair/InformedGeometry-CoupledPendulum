% A Generative Model for Questionnaire-type Data
% as suggested in [J. I. Ankenman, PhD thesis, pp 5--7  Yale University, 2014]
% 
% In this short example, we generate a smooth artificial probability field. We then sample the field as independent
% Bernoulli variables, and then recover the underlying probability field from the sampled data, using binary trees built
% on the known geometry.
close all;
clear all;
clc;
dbstop if error;
n_rows = 100;
n_columns = 100;

means_matrix = zeros(n_rows,n_columns);
Ivec = linspace(0, pi/4, n_rows);
Jvec = linspace(0, 2*pi, n_columns);

[Imat, Jmat] = meshgrid(Ivec, Jvec);
means_matrix = sin((Imat+Jmat+2*Imat.*Jmat)/2.0)*1.0;
rng('shuffle');
pf = means_matrix/2.0+0.5;
orig_data = binornd(1,pf,[n_rows,n_columns]);

params.init_aff.on_rows = true;
params.init_aff.metric = 'cosine_similarity';

params.init_aff.thresh = 0;
params.row_tree.eigs_num = 15;
params.row_tree.constant=2;
params.row_tree.min_joins_percentage = 0.4;
params.col_tree.constant=1;
params.col_tree.min_joins_percentage = 0.5;
params.col_tree.eigs_num = 15;
params.n_iters = 10;
params.col_emd.beta=0;
params.col_emd.alpha = 1;
params.col_emd.eps = 1;
params.row_emd.beta=0;
params.row_emd.alpha = 1;
params.row_emd.eps = 1;

[ row_tree, col_tree ] = RunQuestionnaire( params, orig_data );
order_inds = randperm(n_rows);
 perm_data = orig_data(order_inds,:);
 [ perm_row_tree, perm_col_tree ] = RunQuestionnaire( params, perm_data );
 close all;
reord=[];
parents = unique(perm_row_tree{1}.super_folders);
 for n=1:length(parents)
     inds = find(perm_row_tree{1}.super_folders==parents(n));
     reord = [reord inds];
 end
row_ordered_perm_data = perm_data(reord,:);   

figure;
subplot(2,2,1);
imagesc(means_matrix);
title('Probability Field Means');
subplot(2,2,2);

imagesc(orig_data);
title('Field Realization');



subplot(2,2,3);
imagesc(perm_data);
title('Permuted Realization');
subplot(2,2,4);

imagesc(row_ordered_perm_data);
title('Reordered Permuted Realization');


