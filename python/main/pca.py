import hail as hl
from main.shrink_factor import *


def to_singular_values(eigens):
    return np.sqrt(np.array(eigens)).tolist()


def pca_single(
        mt: hl.MatrixTable,
        n_evecs: int = 10,
        n_evecs_outliers: int = 20,
        maf_thresh: float = 0.05,
        remove_outliers: bool = True,
        sigma_thresh: float = 5.0):

    # count number of samples and number of variants
    ncols = mt.key_cols_by().count_cols()
    nrows = mt.key_cols_by().count_rows()

    # raise error if ncols or nrows are less than number of PCs
    if ncols <= n_evecs or nrows <= n_evecs:
        raise ValueError('Number of samples (%s) or number of variants (%s) are less than %s' % (ncols, nrows, n_evecs))

    # variant QC
    mt = hl.variant_qc(mt)

    # filter by MAF
    mt = mt.filter_rows(
        (mt.variant_qc.AF[0] > maf_thresh) &
        (mt.variant_qc.AF[1] > maf_thresh), keep=True)

    # perform PCA
    eigens, scores_ht, loadings_ht = hl.hwe_normalized_pca(mt.GT, compute_loadings=True, k=n_evecs)

    # compute number of evecs for outlier detection
    n_evecs_outliers = n_evecs if n_evecs <= n_evecs_outliers else n_evecs_outliers

    # calculate the summary statistics for each PC
    evec_mean = []
    evec_stdv = []
    for i in range(n_evecs_outliers):
        stats = scores_ht.aggregate(hl.agg.stats(scores_ht.scores[i]))
        evec_mean.append(stats.mean)
        evec_stdv.append(stats.stdev)

    # filter out outlier samples from matrix table
    if remove_outliers:
        for i in range(n_evecs_outliers):
            scores_ht = scores_ht.filter(
                hl.abs(scores_ht.scores[i] - evec_mean[i]) < sigma_thresh * evec_stdv[i], keep=True)
        mt = mt.filter_cols(hl.is_defined(scores_ht[mt.col_key]), keep=True)

    # convert from eigenvalues to singular values
    singulars = to_singular_values(eigens)

    # calculate m effective and sigma effective
    m_eff, sigma_eff = fit_median_iqr(eigens, m_guess=nrows, sigma_guess=1.0)
    shrink_scores = [shrinkage_scores(eig / sigma_eff, ncols, m_eff) for eig in eigens]
    shrink_scores = [float(s) for s in shrink_scores]

    # annotate loadings
    loadings_ht = (loadings_ht.annotate(AF=mt.rows()[loadings_ht.key].variant_qc.AF)
        .annotate_globals(
        n_evecs=n_evecs,
        evec_mean=evec_mean,
        evec_stdv=evec_stdv,
        sigma_thresh=sigma_thresh,
        eigens=eigens,
        singulars=singulars,
        shrink_scores=shrink_scores))

    # scale scores with singular values
    scores_ht = scores_ht.transmute(scores=[scores_ht.scores[i] * singulars[i] for i in range(n_evecs)])
    scores_ht = (scores_ht
        .annotate_globals(n_evecs=n_evecs
    ))

    # return
    return mt, scores_ht, loadings_ht



def pca(mt: hl.MatrixTable,
        n_evecs: int = 10,
        n_evecs_outliers: int = 20,
        remove_outliers: bool = True,
        n_outlieriters: int = 5,
        sigma_thresh: float = 5.0):

    if remove_outliers:
        for iter in range(n_outlieriters-1):
            mt, scores_ht, loadings_ht = pca_single(
                mt=mt,
                remove_outliers=remove_outliers,
                n_evecs=n_evecs,
                n_evecs_outliers = n_evecs_outliers,
                sigma_thresh=sigma_thresh)

    # do one more round of PCA
    mt, scores_ht, loadings_ht = pca_single(
        mt=mt,
        remove_outliers=remove_outliers,
        n_evecs=n_evecs,
        n_evecs_outliers = n_evecs_outliers,
        sigma_thresh=sigma_thresh)

    # return scores and loadings
    return scores_ht, loadings_ht


def filter_scores_ht(
        scores_ht: hl.Table,
        n_evecs: int = 10,
        n_evecs_location:str= 'n_evecs'):
    n_original_evecs = hl.eval(scores_ht[n_evecs_location])
    if n_original_evecs > n_evecs:
        scores_ht = scores_ht.transmute(scores=scores_ht.scores[0:n_evecs])
        scores_ht = scores_ht.transmute_globals(n_evecs=n_evecs)
    return scores_ht


def filter_loadings_ht(
        loadings_ht: hl.Table,
        n_eves: int = 10,
        n_evecs_location:str= 'n_evecs'):
    n_original_evecs = hl.eval(loadings_ht[n_evecs_location])
    if n_original_evecs > n_eves:
        loadings_ht = loadings_ht.transmute(loadings=loadings_ht.loadings[0:n_eves])
        loadings_ht = loadings_ht.transmute_globals(
            n_evecs=n_eves,
            evec_mean=loadings_ht.evec_mean[0:n_eves],
            evec_stdv=loadings_ht.evec_stdv[0:n_eves],
            eigens=loadings_ht.eigens[0:n_eves],
            singulars=loadings_ht.singulars[0:n_eves],
            shrink_scores=loadings_ht.shrink_scores[0:n_eves])
    return loadings_ht


def pca_project(
        mt:hl.MatrixTable,
        loadings_ht:hl.Table,
        loading_location:str = 'loadings',
        af_location:str = 'AF',
        singulars_location: str = 'singulars',
        shrink_scores_location: str = 'shrink_scores',
        eigenvec_mean_location: str = 'evec_mean',
        eigenvec_stdv_location: str = 'evec_stdv',
        sigma_threshold_location: str = 'sigma_thresh',
        remove_outliers: bool = True,
        correct_shrinkage: bool = False) -> hl.Table:

    singulars = hl.eval(loadings_ht[singulars_location])
    shrink_scores = hl.eval(loadings_ht[shrink_scores_location])
    evec_mean = hl.eval(loadings_ht[eigenvec_mean_location])
    evec_stdv = hl.eval(loadings_ht[eigenvec_stdv_location])
    sigma_thresh = hl.eval(loadings_ht[sigma_threshold_location])

    # count number of eigen-vectors
    n_evecs = len(singulars)
    n_variants = loadings_ht.count()

    # project genotype with the loadings
    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location][1]
    )

    # filter matrix table
    mt = mt.filter_rows(
        hl.is_defined(mt.pca_loadings) &
        hl.is_defined(mt.pca_af) &
        (mt.pca_af > 0) &
        (mt.pca_af < 1))

    # projection
    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))
    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))
    scores_ht = mt.cols().select('scores')

    # remove outliers
    if remove_outliers:
        for i, (m, std) in enumerate(zip(evec_mean, evec_stdv)):
            scores_ht = scores_ht.filter(hl.abs(scores_ht.scores[i] - m) < sigma_thresh * std, keep=True)

    if correct_shrinkage:
        scores_ht = scores_ht.transmute(scores=[scores_ht.scores[i] / shrink_scores[i] * singulars[i] for i in range(n_evecs)])
    else:
        scores_ht = scores_ht.transmute(scores=[scores_ht.scores[i] * singulars[i] for i in range(n_evecs)])

    # return scores
    return scores_ht


def compute_full_spectrum(mt: hl.MatrixTable, n_evecs:int=20,
                          remove_outliers:bool=True, n_outlieriters:int=2, sigma_thresh:float=0.5):

    # count PCA scores and loadings
    scores_ht, loadings_ht = pca(
        mt=mt,
        n_evecs=n_evecs,
        n_evecs_outliers=10,
        remove_outliers=remove_outliers,
        n_outlieriters=n_outlieriters,
        sigma_thresh=sigma_thresh)

    # filter out outlier samples
    mt = mt.filter_cols(
        hl.is_defined(scores_ht[mt.col_key]))
    mt = mt.filter_rows(
        hl.is_defined(loadings_ht[mt.row_key]))

    # HWE-normalize genotypes
    mt = mt.annotate_rows(GT_stats=hl.agg.stats(mt.GT.n_alt_alleles()))
    mt = mt.annotate_rows(hwe_sd=(mt.GT_stats.mean * (1 - mt.GT_stats.mean / 2)) ** 0.5)
    mt = mt.annotate_entries(normed_GT=hl.or_else((mt.GT.n_alt_alleles() - mt.GT_stats.mean) / mt.hwe_sd, 0.0))

    # convert to BlockMatrix
    G = hl.linalg.BlockMatrix.from_entry_expr(mt.normed_GT, block_size=256).to_numpy()

    # number of variants
    n_variants = G.shape[0]

    # Compute eigenvalues & sort in descending
    eigens = np.linalg.eigvalsh((G.T @ G) / n_variants)[::-1]

    return eigens, scores_ht, loadings_ht