function [ params ] = SetQuestParams( row_alpha, row_beta, col_alpha, col_beta )
% Sets parameters for the questionnaire algirithm (user accesible).
% 
% data_file:                   data file path. this file must contain a variable  
%                              named 'data' whose rows correspond to sensors and columns
%                              correspond to observations).
% to_normalize:                boolean variable. if true, data should be normalized.
% normalize_over_rows:         boolean variable. if true, normalization is applied on the rows 
%                              (relavent iff to_normalize is true).
% normalization_type:          choices are: 'by_std' (default)| 'scaling_between_0_and_1'| 'summing_to_unity'
%                              (relavent iff to_normalize is true).
% 
% init_aff_on_rows:            boolean variable. if true, the initial affinity is on the rows. 
% init_aff_metric:             choices are: 'cosine_similarity' (default)| 'gaussian'.
% init_aff_thresh:             conatrains the minimal affinities. positive number, default is 0  
%                              (relavent iff init_aff_type is 'cosine_similarity'). 
% init_aff_eps:                defines locality. positive number, default is 1 
%                              (relavent iff init_aff_type is 'gaussian').
% init_aff_knn:                number of nearest neighbors to consider. positive integer, 
%                              default is 5 (relavent iff init_aff_type is 'gaussian').
% 
% n_iters:                     number of iterations to carry out. positive integer, default is 1. 
% 
% row_constant:                conatrains the size of clusters. positive number, default is 1.
% row_min_joins_percentage:    conatrains the minimal number of joins per level. ranges between 0 and 1. 
% row_alpha:                   real number, default is 1. 
%                              alpha=0: all tree levels are equally weighted.
%                              alpha>0: puts higher weights on levels closer to the root.
%                              alpha<0: puts higher weights on levels closer to the leaves.
% row_beta:                    real number, default is 0. 
%                              beta=0: all tree folders are equally weighted.
%                              beta>0: puts higher weights on larger folders.
%                              beta<0: puts higher weights on smaller folders.
% row_eps:                     defines locality. positive number, default is 1.
% 
% col_constant:                conatrains the size of clusters. positive number, default is 1.
% col_min_joins_percentage:    conatrains the minimal number of joins per level. ranges between 0 and 1. 
% col_alpha:                   real number, default is 1. 
%                              alpha=0: all tree levels are equally weighted.
%                              alpha>0: puts higher weights on levels closer to the root.
%                              alpha<0: puts higher weights on levels closer to the leaves.
% col_beta:                    real number, default is 0. 
%                              beta=0: all tree folders are equally weighted.
%                              beta>0: puts higher weights on larger folders.
%                              beta<0: puts higher weights on smaller folders.
% col_eps:                     defines locality. positive number, default is 1.
% 
% Output: 
%     params - struct array with all user parameters 
%--------------------------------------------------------------------------

data_file = 'data_files/sinusoidal_pf.mat';
to_normalize = false;
normalize_over_rows = true;
normalization_type = 'by_std';

init_aff_on_rows = true;
init_aff_metric = 'cosine_similarity';
init_aff_thresh = 0;
init_aff_eps = 1;
init_aff_knn = 5;

n_iters = 1;

row_constant = 0.5;
row_min_joins_percentage = 0.1;
% row_alpha = 0;
% row_beta = 1;
row_eps = 1;

col_constant = 0.5;
col_min_joins_percentage = 0.1;
% col_alpha = 0;
% col_beta = 1;
col_eps = 1;

%validate arguments values:
% assert(isboolean(to_normalize));
% assert(isboolean(normalize_over_rows));
assert(~to_normalize || strcmp(normalization_type,'by_std') || ...
    strcmp(normalization_type,'scaling_between_0_and_1') || strcmp(normalization_type,'summing_to_unity'));
% assert(isboolean(init_aff_on_rows));
assert(strcmp(init_aff_metric,'cosine_similarity') || strcmp(init_aff_metric,'gaussian'));
assert(~strcmp(init_aff_metric,'cosine_similarity') || (init_aff_thresh >= 0));
assert(~strcmp(init_aff_metric,'gaussian') || (init_aff_eps > 0));
assert(~strcmp(init_aff_metric,'gaussian') || is_natural(init_aff_knn));
assert(is_natural(n_iters));
assert(row_constant > 0);
assert((row_min_joins_percentage >= 0) && (row_min_joins_percentage <= 1));
assert(row_eps > 0);
assert(col_constant > 0);
assert((col_min_joins_percentage >= 0) && (col_min_joins_percentage <= 1));
assert(col_eps > 0);

%--------------------------------------------------------------------------

params.data.file = data_file;
params.data.to_normalize = to_normalize;
params.data.over_rows = normalize_over_rows;
params.data.normalization_type = normalization_type;

params.init_aff.on_rows = init_aff_on_rows;
params.init_aff.metric = init_aff_metric;
params.init_aff.thresh = init_aff_thresh;
params.init_aff.eps = init_aff_eps;
params.init_aff.knn = init_aff_knn;

params.n_iters = n_iters;

params.row_tree.constant = row_constant;
params.row_tree.min_joins_percentage = row_min_joins_percentage;
params.row_emd.alpha = row_alpha;
params.row_emd.beta = row_beta;
params.row_emd.eps = row_eps;

params.col_tree.constant = col_constant;
params.col_tree.min_joins_percentage = col_min_joins_percentage;
params.col_emd.alpha = col_alpha;
params.col_emd.beta = col_beta;
params.col_emd.eps = col_eps;

end


function [ res ] = is_natural( a )
% Determines whether a given number is a natural number.
%--------------------------------------------------------------------------

if (a > 0) && (a == round(a)),
    res = true;
else
    res = false;
end

end


