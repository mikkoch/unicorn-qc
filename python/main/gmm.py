import random
from main.utils import *
from typing import *
from pprint import pprint
from sklearn import mixture
from sklearn.mixture import GaussianMixture

from .utils import to_pandas_and_explode_array

class GMM:
    def __init__(self, n_clusters:int, covariance_type:str= 'full',
                 tol:float=1e-3, reg_covar=1e-6, max_iter:int=100, n_init:int=20, homogeneity_threshold:float=1e-3):
        self.tol = tol
        self.n_init = n_init
        self.max_iter = max_iter
        self.reg_covar = reg_covar
        self.covariance_type=covariance_type
        self.homogeneity_threshold=homogeneity_threshold
        self.n_significant_ancestry = n_clusters - 1
        self.n_clusters = n_clusters


    def fit(self,
            scores: Union[hl.Table, pd.DataFrame],
            pc_cols: Union[hl.expr.ArrayExpression, List[str]]=None,
            output_col:str='cluster',
            #min_prob: float=0.3,
            #missing_label: str='oth',
            out_col_prefix:str='root') -> Union[hl.Table, pd.DataFrame]:

        hail_input = isinstance(scores, hl.Table)

        if hail_input:
            scores_pd = scores.to_pandas()
            scores_pd = expand_pd_array_col(scores_pd, 'scores', out_cols_prefix='PC_{prefix}-'.format(prefix=out_col_prefix))
            pc_cols = [col for col in scores_pd if col.startswith('PC')]
        else:
            scores_pd = scores

        gmm_clf: GaussianMixture = mixture.GaussianMixture(
            n_components=self.n_clusters,
            covariance_type=self.covariance_type,
            tol=self.tol,
            reg_covar=self.reg_covar,
            max_iter=self.max_iter,
            n_init=self.n_init,
            random_state=random.randint(0, 2**32-1))

        gmm_clf = gmm_clf.fit(scores_pd[pc_cols].values)
        scores_pd[output_col] = gmm_clf.predict(scores_pd[pc_cols].values)
        scores_pd[output_col] = out_col_prefix + "-" + scores_pd[output_col].astype(str)
        probs = gmm_clf.predict_proba(scores_pd[pc_cols].values)
        probs = pd.DataFrame(probs, columns=['prob_{prefix}-{p}'.format(p=p, prefix=out_col_prefix) for p in range(self.n_clusters)])
        probs['prob_max-{prefix}'.format(prefix=out_col_prefix)] = probs.max(axis=1)
        scores_pd = pd.concat([scores_pd, probs], axis=1)
        #scores_pd.loc[probs['prob_max-{prefix}'.format(prefix=out_col_prefix)] < min_prob, output_col] = missing_label
        #scores_pd = scores_pd.drop(pc_cols, axis='columns')
        self.gmm_clf = gmm_clf

        if hail_input:
            scores_ht = hl.Table.from_pandas(scores_pd, key=list(scores.key))
            return scores_ht
        else:
            return scores_pd


    def mahalanobis_distance(
            self,
            scores: Union[hl.Table, pd.DataFrame],
            pc_cols: Union[hl.expr.ArrayExpression, List[str]]=None,
            #output_col: str = 'cluster',
            out_col_prefix: str = 'root',
            #max_dist: float = 5,
            #missing_label: str = 'oth',
    ):

        hail_input = isinstance(scores, hl.Table)
        if hail_input:
            scores_pd = scores.to_pandas()
            scores_pd = expand_pd_array_col(scores_pd, 'scores', out_cols_prefix='PC_{prefix}-'.format(prefix=out_col_prefix))
            pc_cols = [col for col in scores_pd if col.startswith('PC')]
        else:
            scores_pd = scores

        means = self.gmm_clf.means_
        precisions = self.gmm_clf.precisions_
        n_components, _ = means.shape
        n_samples, _ = scores_pd[pc_cols].shape

        distances = np.empty((n_samples, n_components))

        for i in range(n_samples):
            for k in range(n_components):
                x = scores_pd[pc_cols].values[i, :]
                diff = x - means[k]
                distances[i, k] = np.sqrt(np.dot(np.matmul(diff.transpose(), precisions[k]), diff))

        distances = pd.DataFrame(
            data=distances,
            columns=['dist_{prefix}-{p}'.format(p=p, prefix=out_col_prefix) for p in range(self.n_clusters)])
        distances['dist_min-{prefix}'.format(prefix=out_col_prefix)] = distances.min(axis=1)
        scores_pd = pd.concat([scores_pd, distances], axis=1)
        #scores_pd.loc[distances['dist_min-{prefix}'.format(prefix=out_col_prefix)] > max_dist, output_col] = missing_label
        #scores_pd = scores_pd.drop(pc_cols, axis='columns')

        if hail_input:
            scores_ht = hl.Table.from_pandas(scores_pd, key=list(scores.key))
            return scores_ht
        else:
            return scores_pd



    def predict(self,
                scores:hl.Table,
                pc_cols: Union[hl.expr.ArrayExpression, List[str]] = None,
                output_col: str = 'cluster',
                #min_prob: float = 0.3,
                #missing_label: str = 'oth',
                out_col_prefix: str = 'root')->Union[hl.Table, pd.DataFrame]:

        hail_input = isinstance(scores, hl.Table)
        if hail_input:
            scores_pd = scores.to_pandas()
            scores_pd = expand_pd_array_col(scores_pd, 'scores', out_cols_prefix='PC_{prefix}-'.format(prefix=out_col_prefix))
            pc_cols = [col for col in scores_pd if col.startswith('PC')]
        else:
            scores_pd = scores

        scores_pd[output_col] = self.gmm_clf.predict(scores_pd[pc_cols].values)
        scores_pd[output_col] = out_col_prefix + "-" + scores_pd[output_col].astype(str)
        probs = self.gmm_clf.predict_proba(scores_pd[pc_cols].values)
        probs = pd.DataFrame(probs, columns=['prob_{prefix}-{p}'.format(p=p, prefix=out_col_prefix) for p in
                                             range(self.n_clusters)])
        probs['prob_max-{prefix}'.format(prefix=out_col_prefix)] = probs.max(axis=1)
        scores_pd = pd.concat([scores_pd, probs], axis=1)
        #scores_pd.loc[probs['prob_max-{prefix}'.format(prefix=out_col_prefix)] < min_prob, output_col] = missing_label
        #scores_pd = scores_pd.drop(pc_cols, axis='columns')

        if hail_input:
            scores_ht = hl.Table.from_pandas(scores_pd, key=list(scores.key))
            return scores_ht
        else:
            return scores_pd



