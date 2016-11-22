% In this example we generate a smooth artificial probability field.
% then we sample it as independent Bernoulli variables, and attempt to  
% recover the underlying field from the measurements, using the Questionnaire 
% algorithm
% =========================================================================
close all; clear; clc;


%% ========================================================================
% create input data:
% ========================================================================
n_rows = 300;
n_cols = 150;

[x, y] = meshgrid((1:n_cols)*2*pi/n_cols, (1:n_rows)*2*pi/n_rows);
pf = 0.5*(1+sin((x+y+2*x.*y)/2));
figure, imagesc(pf), colormap jet, axis on, title('Probability Field'), colorbar

orig_data = binornd(1,pf); 
orig_data(orig_data == 0) = -1;
figure, imagesc(orig_data), colormap gray, axis on, title('Field Realization'), colorbar

row_perm = randperm(n_rows);
col_perm = randperm(n_cols);
data = orig_data(row_perm,:); 
data = data(:,col_perm);
figure, imagesc(data), colormap gray, axis on, title('Shuffled Data'), colorbar


%% ========================================================================
% set parameters and run Questionnaire:
% ========================================================================
row_alpha = 0;      
col_alpha = 0;      
row_beta = 1;
col_beta = 1;

params = SetQuestParams(row_alpha, row_beta, col_alpha, col_beta);
[row_tree, col_tree] = RunQuestionnaire(params, data);


%% ========================================================================
% rearrange the data based on the induced metrics:
% ========================================================================
row_aff = CalcEmdAff(data.', col_tree, params.row_emd);
col_aff = CalcEmdAff(data, row_tree, params.col_emd);

row_thresh = 0.6; 
col_thresh = 0.6;
% fairly high thresholds on the affinities (data-dependent), which are aimed  
% to collapse the diffusion map into a curve.
row_aff = threshold(row_aff, row_thresh);
col_aff = threshold(col_aff, col_thresh);
OrganizeData(orig_data, data, row_aff, col_aff, row_perm, col_perm);




