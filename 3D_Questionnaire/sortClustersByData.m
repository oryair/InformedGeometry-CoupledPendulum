function clusteringRes = sortClustersByData(data, clustering)
clusters = unique(clustering);

L = length(unique(clustering));
if L < 9
    ordering = perms(1:L);
    for o = 1:size(ordering, 1)
        clustering1 = zeros(size(clustering));
        for ci = 1:length(clusters)
            clustering1(clustering==  clusters(ci)) = ordering(o, ci);
        end
        [~, orderingAll] = sort(clustering1, 'ascend');
        
        err(o) = sum(sum(abs(data(orderingAll, :)-data)));
    end
    [~,m]=min(err);
    clusteringRes = zeros(size(clustering));
    for ci = 1:length(clusters)
        clusteringRes(clustering==  clusters(ci)) = ordering(m, ci);
    end
    
else
    [~, ordering1] = sort(clustering, 'descend');
    [~, ordering2] = sort(clustering, 'ascend');
    clustersRev = flipud(unique(clusters));
    clusteringRes = -clustering;
    if sum(sum(abs(data(ordering1, :)-data))) < sum(sum(abs(data(ordering2, :)-data)))
        for ci = 1:length(clusters)
            clusteringRes(clustering==  clusters(ci)) = clustersRev(ci);
        end
    else
        clusteringRes = -clusteringRes;
    end
end