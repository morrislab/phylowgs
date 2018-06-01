function ClusterUtil() {
}

ClusterUtil.separate_clusters_by_size = function(clusters, cluster_size_criteria){
    if(Config.report_small_clusters){
        return {small: [], large: clusters};
    }
    var small_clusters = {};
    var large_clusters = {};
    Object.keys(clusters).forEach(function(cidx){
        var C = clusters[cidx];
        if(C.members.length < cluster_size_criteria){
        small_clusters[cidx] = C;
        }else{
        large_clusters[cidx] = C;
        }
    })
    return {small: small_clusters, large: large_clusters};
}
