from main.n_spikes import *
from main.tree import *

def build_tree(
        mt_path: str,
        loadings_ht: hl.Table = None,
        scores_ht: hl.Table = None,
        tree_meta: Tree = None,
        tree_branch: str = "root",
        remove_outliers: bool = True,
        n_outlieriters: int = 2,
        sigma_thresh: float = 5.0,
        create_binary_tree: bool = False,
        n_eves_for_regression:int=100,
        n_samples_thresh:int=200,
        studentized_residual_threshold:float=1e-5) -> (Tree, hl.Table, hl.Table):

    mt = hl.read_matrix_table(mt_path)

    # create Tables for storing returned values, and select samples for tree building
    if loadings_ht is None: loadings_ht = mt.rows().select()
    if tree_meta is None: tree_meta = Tree()
    if scores_ht is None:
        scores_ht = mt.cols().select()
    else:
        cluster_location = 'cluster-{branch}'.format(
            branch = Node.get_parent(tree_branch))
        mt = mt.filter_cols(
            (hl.is_defined(scores_ht[mt.col_key][cluster_location])) &
            (scores_ht[mt.col_key][cluster_location] == tree_branch))

    print("Use {n_samples} samples and {n_variants} variants for building ancestry tree at branch {branch}".format(
        n_samples=mt.count_cols(),
        n_variants=mt.count_rows(),
        branch=tree_branch))

    n_significant_eigens, scores, loadings = compute_n_significant_eigens(
        mt=mt,
        remove_outliers=remove_outliers,
        n_outlieriters=n_outlieriters,
        sigma_thresh=sigma_thresh,
        n_eves_for_regression=n_eves_for_regression,
        sample_size_threshold=n_samples_thresh,
        studentized_residual_threshold=studentized_residual_threshold)
    is_leaf = n_significant_eigens == 0
    if create_binary_tree and n_significant_eigens > 0: n_significant_eigens = 1
    print("There are {n} significant eigenvalues at branch {branch}".format(
        n=n_significant_eigens,
        branch=tree_branch))

    if is_leaf:
        cluster_meta = Node(
            name=tree_branch,
            n_populations=0,
            n_samples=mt.count_cols(),
            samples=mt.s.collect(),
            is_leaf=is_leaf)
        tree_meta.add_node(cluster_meta)
        return tree_meta, scores_ht, loadings_ht

    else:
        # Filter to spikes
        scores = filter_scores_ht(
            scores_ht=scores,
            n_evecs=n_significant_eigens)
        loadings = filter_loadings_ht(
            loadings_ht=loadings,
            n_eves=n_significant_eigens)

        # Perform GMM clustering
        gmm_clf = GMM(n_clusters=n_significant_eigens+1)
        cluster_ht = gmm_clf.fit(
            scores=scores,
            output_col='cluster-{postfix}'.format(postfix=tree_branch),
            out_col_prefix=tree_branch)
        cluster_ht = gmm_clf.mahalanobis_distance(
            scores=cluster_ht,
            out_col_prefix=tree_branch)

        # Create node meta information
        cluster_meta = Node(
            name=tree_branch,
            n_populations=n_significant_eigens+1,
            n_samples=mt.count_cols(),
            samples=cluster_ht.s.collect(),
            gmm=gmm_clf,
            is_leaf=is_leaf)
        tree_meta.add_node(cluster_meta)

        # Annotate results
        scores_ht = scores_ht.annotate(
            **cluster_ht[scores_ht.key])
        scores_ht = scores_ht.annotate(
            **{'scores_{branch}'.format(branch=tree_branch): cluster_ht[scores_ht.key].scores})

        loadings_ht = loadings_ht.annotate(**{
            'loadings_{branch}'.format(branch=tree_branch): loadings[loadings_ht.key].loadings,
            'AF_{branch}'.format(branch=tree_branch): loadings[loadings_ht.key].AF})
        loadings_ht = loadings_ht.annotate_globals(**{
            'n_evecs-{branch}'.format(branch=tree_branch): hl.eval(loadings.n_evecs),
            'evec_mean-{branch}'.format(branch=tree_branch): hl.eval(loadings.evec_mean),
            'evec_stdv-{branch}'.format(branch=tree_branch): hl.eval(loadings.evec_stdv),
            'sigma_thresh-{branch}'.format(branch=tree_branch): hl.eval(loadings.sigma_thresh),
            'eigens-{branch}'.format(branch=tree_branch): hl.eval(loadings.eigens),
            'singulars-{branch}'.format(branch=tree_branch): hl.eval(loadings.singulars),
            'shrink_scores-{branch}'.format(branch=tree_branch): hl.eval(loadings.shrink_scores)})

        # Iterate
        n_clusters = gmm_clf.gmm_clf.n_components
        for idx_cluster in range(n_clusters):
            child_branch = '{parent_branch}-{child_index}'.format(
                parent_branch=tree_branch,
                child_index=idx_cluster)

            tree_meta, scores_ht, loadings_ht = build_tree(
                mt_path=mt_path,
                scores_ht=scores_ht,
                loadings_ht=loadings_ht,
                tree_meta=tree_meta,
                tree_branch=child_branch,
                remove_outliers=remove_outliers,
                n_outlieriters=n_outlieriters,
                n_eves_for_regression=n_eves_for_regression,
                n_samples_thresh=n_samples_thresh,
                studentized_residual_threshold=studentized_residual_threshold,
                sigma_thresh=sigma_thresh,
                create_binary_tree=create_binary_tree)

            scores_ht = scores_ht.cache()
            loadings_ht = loadings_ht.cache()

        return tree_meta, scores_ht, loadings_ht
