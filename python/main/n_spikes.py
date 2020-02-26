from main.pca import *
from pprint import pprint
from typing import *
import statsmodels.api as sm
from scipy.stats import t


def n_spikes(
        s: Union[List[float], np.ndarray],
        n_evecs_for_reg: int,
        n_samples:int,
        sample_size_threshold:int=100,
        cook_distance_threshold: float=1.0,
        p_threshold:float=1e-5,
        method='residual'):

    if isinstance(s, np.ndarray):
        s = list(s)

    if n_samples < sample_size_threshold:
        return 0

    else:
        for i in range(20):
            Y = np.array(s[i: (i + n_evecs_for_reg)])
            X = np.array([((index + 0.5)/(n_samples-i)) ** (2.0/3.0) for index in range(n_evecs_for_reg)])
            lm = sm.OLS(Y, sm.add_constant(X)).fit()

            # Calculate studentized residuals and cook's distance
            influence = lm.get_influence()
            cooks_d = influence.cooks_distance
            studentized_residuals = influence.resid_studentized_external

            # Threshold
            if method == 'residual':
                p = 1-t.cdf(studentized_residuals, df=sample_size_threshold-3, loc=0, scale=1)
                if p[0] >= p_threshold: break
            elif method == 'cook':
                if cooks_d[0][0] < cook_distance_threshold: break
            else:
                raise ValueError(
                    'Invalid method input: {method}'.format(method=method))
        return i


def n_spikes_studentized(s, n_samples, n=100, p_threshold=1e-5):
    for i in range(20):
        Y = s[i:i+n]
        X = ((np.arange(n) + 0.5) / (n_samples - i)) ** (2 / 3)
        lm = sm.OLS(Y, sm.add_constant(X)).fit()

        # Calculate studentized residuals
        influence = lm.get_influence()
        studentized_res = influence.resid_studentized_external
        if studentized_res[0] < 0: break

        p = 1 - t.cdf(studentized_res, df=n-3, loc=0, scale=1)
        if p[0] >= p_threshold: break

    return i


def compute_n_significant_eigens(
        mt: hl.matrixtable,
        n_evecs:int=20,
        remove_outliers:bool=True,
        n_outlieriters:int=2,
        sigma_thresh:float=5.0,
        n_eves_for_regression:int=100,
        sample_size_threshold:int=200,
        studentized_residual_threshold:float=1e-5):

    n_samples = mt.count_cols()

    if n_samples < sample_size_threshold:
        return 0, None, None
    else:
        s, scores_ht, loadings_ht = compute_full_spectrum(
            mt=mt,
            n_evecs=n_evecs,
            remove_outliers=remove_outliers,
            n_outlieriters=n_outlieriters,
            sigma_thresh=sigma_thresh)

        n_samples = scores_ht.count()
        n_significant_eigens = n_spikes_studentized(
            s=s,
            n_samples=n_samples,
            n=n_eves_for_regression,
            p_threshold=studentized_residual_threshold)

        # print #significant eigen-vectors
        pprint("Number of significant eigen-vectors: %d" %n_significant_eigens)

        # return
        return n_significant_eigens, scores_ht, loadings_ht
