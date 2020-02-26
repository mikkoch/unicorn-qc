import hail as hl
from main.tree import *

def fine_matching(
        cluster_co_ht: hl.Table,
        cluster_ca_ht: hl.Table,
        scores_co_ht: hl.Table = None,
        scores_ca_ht: hl.Table = None,
        sample_size_threshold: int = 0,
        cluster_col: str = 'cluster'):

    cluster_co_ht = cluster_co_ht.annotate(is_case=0)
    cluster_ca_ht = cluster_ca_ht.annotate(is_case=1)
    cluster_ht = cluster_co_ht.union(cluster_ca_ht)

    # Group by Case Control Status
    cluster_summ_ht = (cluster_ht
        .group_by(cluster_ht[cluster_col])
        .aggregate(n_cases=hl.agg.sum(cluster_ht.is_case == 1),
                   n_conts=hl.agg.sum(cluster_ht.is_case == 0)))

    # filter by sample size threshold
    cluster_summ_ht = cluster_summ_ht.filter(
        (cluster_summ_ht.n_cases > sample_size_threshold) &
        (cluster_summ_ht.n_conts > sample_size_threshold),
        keep=True)

    # cluster labels
    clusters = hl.array(cluster_summ_ht.cluster.collect())

    # filter cluster Table
    cluster_ht = cluster_ht.filter(
        clusters.contains(cluster_ht.cluster), keep=True)

    # throw away non-PCs columns
    if scores_ca_ht is not None:
        cols = list(scores_ca_ht.row_value.keys())
        pc_cols = [col for col in cols if col.startswith('PC')]
        scores_ca_ht = scores_ca_ht.select(*pc_cols)
    if scores_co_ht is not None:
        cols = list(scores_co_ht.row_value.keys())
        pc_cols = [col for col in cols if col.startswith('PC')]
        scores_co_ht = scores_co_ht.select(*pc_cols)

    #
    # for cluster in clusters:
    #     print(cluster)

    scores_ht = scores_co_ht.union(scores_ca_ht)

    scores_ht = scores_ht.annotate(
        cluster=cluster_ht[scores_ht.key].cluster,
        is_case=cluster_ht[scores_ht.key].is_case)
    scores_ht = scores_ht.filter(
        hl.is_missing(scores_ht.cluster), keep=False)


    return cluster_ht, scores_ht