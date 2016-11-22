function params = SetGenericDimsQuestParams(dims, isTopDown)


VERBOSE_NO_OUTPUT      = 0; %#ok
VERBOSE_PROPMT_OUTPUT  = 1; %#ok
VERBOSE_FIGURES_OUTPUT = 2; %#ok

AFF_EUCLID_METRIC               = 'euc';                       %#ok
AFF_COSINE_SIMILARITY           = 'cosine_similarity';         %#ok
AFF_COSINE_SIMILARITY_ON_TRAILS = 'cosine_similarityOnTrials'; %#ok


params.n_iters           = 1;
if ~isunix
    params.verbose           = VERBOSE_FIGURES_OUTPUT; % 2 - for printing
else
    params.verbose       = VERBOSE_PROPMT_OUTPUT; % 2 - for printing
end
params.data.to_normalize = false;

for dim_i = 1:dims
    params.tree{dim_i}.runOnEmbdding = true; % build tree based on distances between embeddings (true) or on the given matrix (false)
    params.tree{dim_i}.eigs_num      = 10;
    
    if ~isunix
        params.tree{dim_i}.verbose           = VERBOSE_FIGURES_OUTPUT; % 2 - for printing
    else
        params.tree{dim_i}.verbose       = VERBOSE_PROPMT_OUTPUT; % 2 - for printing
    end
    params.init_aff{dim_i}.metric = AFF_EUCLID_METRIC;
    params.init_aff{dim_i}.knn    = 5;
    params.init_aff{dim_i}.eps    = 1;
    params.init_aff{dim_i}.thresh = 0;
    
    
    
    params.init_aff{dim_i}.initAffineFun = @CalcInitAff;
    params.tree{dim_i}.CalcAffFun        = @CalcEmdAff;
        
    params.emd{dim_i}.beta  = 1;
    params.emd{dim_i}.alpha = 0;
    params.emd{dim_i}.eps   = 1;
    
    if isTopDown
        % relevant only for TD trees:
        params.tree{dim_i}.treeDepth      = 4; % limiting the tree's groth
        params.tree{dim_i}.clusteringAlgo = @svdClassWrapper;
        params.tree{dim_i}.splitsNum      = 9; % how many splits we want at each node
        params.tree{dim_i}.min_cluster    = 12;
        params.tree{dim_i}.buildTreeFun   = @BuildGenericTdTreesViaClustering;
    else
        params.tree{dim_i}.buildTreeFun         = @BuildFlexTree;
        params.tree{dim_i}.runOnEmbdding        = false; % for this option we want false because BuildFlexTree embedd the data
        params.tree{dim_i}.threshold            = 0;
        params.tree{dim_i}.k                    = 2;
        params.tree{dim_i}.min_cluster          = 12;
        params.tree{dim_i}.constant             = .5;
        params.tree{dim_i}.min_joins_percentage = 0.1;
    end
    
end
end





