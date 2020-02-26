from main.tree import *


def get_ca_probs(
        pops_mt: hl.MatrixTable,
        scores_ht:hl.Table,
        rows_cluster_name:str='cluster',
        distance_thresh: float = 100,
        posterior_thresh: float = 0):

    scores_ht_cols = list(scores_ht.row_value.keys())
    posterior_cols = [col for col in scores_ht_cols if col.startswith('prob_max')]
    distance_cols = [col for col in scores_ht_cols if col.startswith('dist_min')]

    scores_ht = scores_ht.annotate(is_outlier=False)

    for distance_col in distance_cols:
        scores_ht = scores_ht.transmute(
            is_outlier=hl.cond(
                ((scores_ht[distance_col] > distance_thresh) | (scores_ht.is_outlier)), True, False))

    for posterior_col in posterior_cols:
        scores_ht = scores_ht.transmute(
            is_outlier=hl.cond(
                ((scores_ht[posterior_col] < posterior_thresh) | (scores_ht.is_outlier)), True, False))

    print("Outlier distribution: {}".format(
        ", ".join(
            f'{pop}: {count}' for pop, count in scores_ht.aggregate(hl.agg.counter(scores_ht.is_outlier)).items()
        )))

    # filter outliers
    scores_ht = scores_ht.filter(scores_ht.is_outlier == True, keep=False)

    # returned Table
    probs_ht = scores_ht.select()
    population_cols = pops_mt[rows_cluster_name].collect()

    for population in population_cols:
        probs_ht = probs_ht.annotate(**{
            population: hl.cond(
                hl.is_defined(
                    scores_ht[probs_ht.key]['prob_{population}'.format(population=population)]),
                    scores_ht[probs_ht.key]['prob_{population}'.format(population=population)],
                    0)})

        parent_pop_iter = Node.get_parent(population)
        while not Node.is_root(parent_pop_iter):
            probs_ht = probs_ht.transmute(**{
                population: hl.cond(
                    hl.is_defined(
                        scores_ht[probs_ht.key]['prob_{population}'.format(population=parent_pop_iter)]),
                    probs_ht[population] * scores_ht[probs_ht.key]['prob_{population}'.format(population=parent_pop_iter)],
                        0)})
            parent_pop_iter = Node.get_parent(parent_pop_iter)

    return probs_ht


def get_cluster_assignment_by_posterior_probability(
        pops_mt:hl.MatrixTable,
        scores_ht:hl.Table,
        rows_cluster_name:str='cluster',
        output_rows_cluster_name:str='cluster',
        distance_thresh: float=100,
        posterior_thresh: float=0):

    population_cols = pops_mt[rows_cluster_name].collect()

    probs_ht = get_ca_probs(
        pops_mt=pops_mt,
        scores_ht=scores_ht,
        rows_cluster_name=rows_cluster_name,
        posterior_thresh=posterior_thresh,
        distance_thresh=distance_thresh)

    probs_ht = probs_ht.annotate(
        probs=[probs_ht[pop_col] for pop_col in population_cols])
    probs_ht = probs_ht.annotate(
        probs_max=hl.max(probs_ht.probs),
        probs_argmax=hl.argmax(probs_ht.probs))

    cluster_pd = probs_ht.select(probs_ht.probs_argmax).to_pandas()
    cluster_pd[output_rows_cluster_name] = pd.Series(
        [population_cols[idx] for idx in cluster_pd.probs_argmax.tolist()])
    cluster_pd = cluster_pd.drop('probs_argmax', axis='columns')
    cluster_ht = hl.Table.from_pandas(cluster_pd, key=list(probs_ht.key))

    return cluster_ht


def get_cluster_assignment_by_tree_meta(
        scores_ht: hl.Table,
        tree_meta: Tree,
        distance_thresh: float = 100,
        posterior_thresh: float = 0,
        n_samples_thresh:int = 100,
        output_rows_cluster_name:str='cluster'):

    scores_ht_cols = list(scores_ht.row_value.keys())
    posterior_cols  = [col for col in scores_ht_cols if col.startswith('prob_max')]
    distance_cols = [col for col in scores_ht_cols if col.startswith('dist_min')]

    scores_ht = scores_ht.annotate(is_outlier=False)

    for distance_col in distance_cols:
        scores_ht = scores_ht.transmute(
            is_outlier=hl.cond(
                (((scores_ht[distance_col] > distance_thresh) & hl.is_defined(scores_ht[distance_col])) | (scores_ht.is_outlier) ), True, False))

    for posterior_col in posterior_cols:
        scores_ht = scores_ht.transmute(
            is_outlier=hl.cond(
                (((scores_ht[posterior_col] < posterior_thresh) & hl.is_defined(scores_ht[posterior_col])) | (scores_ht.is_outlier)), True, False))

    print("Outlier distribution: {}".format(
        ", ".join(
            f'{pop}: {count}' for pop, count in scores_ht.aggregate(hl.agg.counter(scores_ht.is_outlier)).items()
            )))

    # filter outliers
    scores_ht = scores_ht.filter(scores_ht.is_outlier == True, keep=False)

    cluster_ht = tree_meta.get_assignment(
        output_col=output_rows_cluster_name,
        n_samples_thresh=n_samples_thresh,
        hail_output=True)

    cluster_ht = cluster_ht.filter(
        hl.is_defined(scores_ht[cluster_ht.key]), keep=True)

    return cluster_ht

