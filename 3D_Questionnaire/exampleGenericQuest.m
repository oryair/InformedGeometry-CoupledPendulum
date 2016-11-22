%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A simple example of usage of the quest.
% Running on rundom data - just to explain the basic confuguration and
% profiling.
% 
% Written by: Hadas Benisty, 19/7/2016
% hadas.benisty@gmail.com
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2D Small Data
X2d = rand(50);
% get the default for top down trees using pca clustering
params = SetGenericDimsQuestParams(ndims(X2d), true);
% setting the splits to 2 to get a binary tree
for dim_i = 1:ndims(X2d)
    params.tree{dim_i}.splitsNum = 2;
end
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire(params, X2d);
% bottom up
params = SetGenericDimsQuestParams(ndims(X2d), false);
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire(params, X2d);

%% 3D Data
X3d = rand(50, 50, 50);
% get the default for top down trees using pca clustering
params = SetGenericDimsQuestParams(ndims(X3d), true);
% setting the splits to 2 to get a binary tree
for dim_i = 1:ndims(X3d)
    params.tree{dim_i}.splitsNum = 2;
end
% bottom up
params = SetGenericDimsQuestParams(ndims(X3d), false);
params.verbose = 2; % see all figures created by the q.
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire(params, X3d);

%% 2D Big Data 
X2dbig = rand(2500, 50);
% get the default for top down trees using pca clustering
params = SetGenericDimsQuestParams(ndims(X2dbig), true);

for dim_i = 1:ndims(X2dbig)
    params.tree{dim_i}.splitsNum = 2;% setting the splits to 2 to get a binary tree
    params.tree{dim_i}.runOnEmbdding = false; % run on the dist and not on embedding
end
params.verbose = 2; % see all figures created by the q.
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire(params, X2dbig);
% bottom up
params = SetGenericDimsQuestParams(ndims(X2dbig), false);
[ Trees, dual_aff, init_aff, embedding ] = RunGenericDimsQuestionnaire(params, X2dbig);
