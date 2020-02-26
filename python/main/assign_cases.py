from main.pca import *
from main.tree import *

def assign_cases(
        mt_path: str,
        loadings_ht: hl.Table,
        tree_meta: Tree,
        scores_ht: hl.Table = None,
        tree_branch: str = "root",
        remove_outliers:bool = True,
        correct_shrinkage: bool = False,
        min_prob: float = 0,
        max_dist: float = 5.0):

    mt = hl.read_matrix_table(mt_path)
    is_leaf = tree_meta.get_node(tree_branch).is_leaf
    gmm_clf = tree_meta.get_node(tree_branch).get_classifier()

    # create Tables for storing returned values, and select samples for tree building
    if scores_ht is None:
        scores_ht = mt.cols().select()

    print("Assign {n_samples} samples and {n_variants} at branch {branch}".format(
        n_samples=mt.count_cols(),
        n_variants=mt.count_rows(),
        branch=tree_branch))

    if is_leaf:
        return scores_ht

    else:
        # pca project
        scores = pca_project(
            mt=mt,
            loadings_ht=loadings_ht,
            remove_outliers=remove_outliers,
            correct_shrinkage=correct_shrinkage,
            loading_location='loadings_{branch}'.format(branch=tree_branch),
            af_location='AF_{branch}'.format(branch=tree_branch),
            singulars_location='singulars-{branch}'.format(branch=tree_branch),
            shrink_scores_location='shrink_scores-{branch}'.format(branch=tree_branch),
            eigenvec_mean_location='evec_mean-{branch}'.format(branch=tree_branch),
            eigenvec_stdv_location='evec_stdv-{branch}'.format(branch=tree_branch),
            sigma_threshold_location='sigma_thresh-{branch}'.format(branch=tree_branch))

        # Classify
        cluster_ht = gmm_clf.predict(
            scores=scores,
            output_col='cluster-{postfix}'.format(postfix=tree_branch),
            out_col_prefix=tree_branch)
        cluster_ht = gmm_clf.mahalanobis_distance(
            scores=cluster_ht,
            out_col_prefix=tree_branch)

        # Annotate results
        scores_ht = scores_ht.annotate(
            **cluster_ht[scores_ht.key])
        scores_ht = scores_ht.filter(
            hl.is_missing(scores_ht['cluster-{branch}'.format(branch=tree_branch)]), keep=False)

        # Iterate
        n_clusters = gmm_clf.gmm_clf.n_components
        for idx_cluster in range(n_clusters):
            child_branch = '{parent_branch}-{child_index}'.format(
                parent_branch=tree_branch,
                child_index=idx_cluster)

            scores_ht = assign_cases(
                mt_path=mt_path,
                scores_ht=scores_ht,
                loadings_ht=loadings_ht,
                tree_meta=tree_meta,
                tree_branch=child_branch,
                remove_outliers=remove_outliers,
                correct_shrinkage=correct_shrinkage,
                min_prob=min_prob,
                max_dist=max_dist)

            scores_ht = scores_ht.cache()

    return scores_ht




