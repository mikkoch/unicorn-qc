from pprint import pprint
import hail as hl
from collections import *
from typing import *
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import random


def calculate_per_chr_scrt(mt: hl.MatrixTable, scrt_col:str= 'scrt_per_chr')->hl.MatrixTable:
    """
        This function calculate per-chromosome sample call rate, minimum sample call rate cross chromosome
            and annotate original matrix table
        :param MatrixTable mt: Input Hail MatrixTable
        :param str scrt_col: Columns storing minimum sample call rate among chromosome
        :rtype: MatrixTable
        """
    mt = mt.annotate_rows(contig=mt.locus.contig)
    mt_per_chr = mt.group_rows_by(mt.contig).aggregate(call_rate=hl.agg.fraction(hl.is_defined(mt.GT)))
    mt_per_chr = mt_per_chr.annotate_cols(
        call_rate_per_chr=hl.agg.min(mt_per_chr.call_rate))
    mt = mt.annotate_cols(**{scrt_col: mt_per_chr.cols()[mt.col_key].call_rate_per_chr})
    return(mt)


def calculate_inbreeding_coefficients(mt:hl.MatrixTable, ib_col:str= "ib", ib_global:str= "ib_stats", mafrsh:float=1e-8)->hl.MatrixTable:
    """
        This function calculate inbreeding coefficient for each sample, and its summary statistics
        :param MatrixTable mt: Input Hail MatrixTable
        :param str ib_col: Columns storing Inbreeding Coefficients
        :param str ib_global: Global annotation storing summary statistics for F-statistics of IB
        :param float mafrsh: MAF threshold for selecting variants used for calculating inbreeding coefficients
        :rtype: MatrixTable
        """
    mt_cleaned = hl.variant_qc(mt)

    # Filter high quality variants for calculating inbreeding coefficient
    mt_cleaned = mt_cleaned.filter_rows(
        (mt_cleaned.variant_qc.AC[0] > 1) &
        (mt_cleaned.variant_qc.AC[1] > 1) &
        (mt_cleaned.variant_qc.AF[0] > mafrsh) &
        (mt_cleaned.variant_qc.AF[1] > mafrsh) &
        (mt_cleaned.locus.in_autosome()), keep=True)

    # Calculate Inbreeding Coefficients
    mt_cleaned = mt_cleaned.annotate_cols(IB=hl.agg.inbreeding(mt_cleaned.GT, mt_cleaned.variant_qc.AF[1]))
    inbreeding_coef_stats = mt_cleaned.aggregate_cols(hl.agg.stats(mt_cleaned.IB.f_stat))
    mt = mt.annotate_cols(**{ib_col: mt_cleaned.cols()[mt.col_key].IB})
    mt = mt.annotate_globals(**{ib_global: inbreeding_coef_stats})
    return mt


def identify_related_samples(mt:hl.MatrixTable, pihat_threshold:float, ld_pruned:bool=True)->hl.Table:
    """
            This function calculate inbreeding coefficient for each sample, and its summary statistics
            :param MatrixTable mt: Input Hail MatrixTable
            :param str ib_col: Columns storing Inbreeding Coefficients
            :param str ib_global: Global annotation storing summary statistics for F-statistics of IB
            :param float mafrsh: MAF threshold for selecting variants used for calculating inbreeding coefficients
            :rtype: MatrixTable
            """
    if not ld_pruned:
        mt = ld_prune(mt=mt)
    pairs = hl.identity_by_descent(mt, min=pihat_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
    return related_samples_to_remove


def ld_prune(mt:hl.MatrixTable, crsh: float=0.999, mafrsh: float=0.05, hwersh:float=1e-4, r2:float=0.2, pruned_variants_list:bool=False):

    # Report number of samples and number of variants before cleaning
    ncols = mt.count_cols()
    nrows = mt.count_rows()
    print("#variants = %d, #samples = %d before LD pruning" % (nrows, ncols))

    # Variant QC
    mt = hl.variant_qc(mt)

    # Drop all of the SNPs in the region of the LCT locus (hg19 coordinates)
    mt = mt.filter_rows(
        (mt.locus.contig == '2') &
        (mt.locus.position > 129883530) &
        (mt.locus.position < 140283530), keep=False)
    # Drop all of the SNPs in the major histocompatibility complex (hg19 coordinates)
    mt = mt.filter_rows(
        (mt.locus.contig == '6') &
        (mt.locus.position > 24092021) &
        (mt.locus.position < 38892022), keep=False)
    # Drop all of the SNPs in the inverted regions on chromosomes 8 and 17
    mt = mt.filter_rows(
        (mt.locus.contig == '8') &
        (mt.locus.position > 6612592) &
        (mt.locus.position <13455629), keep=False)
    mt = mt.filter_rows(
        (mt.locus.contig == '17') &
        (mt.locus.position > 40546474) &
        (mt.locus.position < 44644684), keep=False)

    # Drop all non-autosomal SNPs
    mt = mt.filter_rows(mt.locus.in_autosome(), keep=True)

    # Filter by minor allele frequencies & call rate
    mt = mt.filter_rows(
        (mt.variant_qc.call_rate > crsh) &
        (mt.variant_qc.AF[0] > mafrsh) &
        (mt.variant_qc.AF[1] > mafrsh) &
        (mt.variant_qc.p_value_hwe > hwersh), keep=True)


    pruned_variant_table = hl.ld_prune(mt.GT, r2=r2)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

    if pruned_variants_list:
        print("{n_variants} number of variants left after LD-pruning".format(
            n_variants=pruned_variant_table.count()))
        return(pruned_variant_table)
    else:
        print("#variants = {n_variants}, #samples = {n_samples} after LD pruning".format(
            n_variants=pruned_mt.count_rows(), n_samples=pruned_mt.count_cols()))
        return(pruned_mt)


def expand_pd_array_col(
        df: pd.DataFrame,
        array_col: str,
        num_out_cols: int = 0,
        out_cols_prefix=None,
        out_1based_indexing: bool = True
) -> pd.DataFrame:
    """
    Expands a Dataframe column containing an array into multiple columns.
    :param DataFrame df: input dataframe
    :param str array_col: Column containing the array
    :param int num_out_cols: Number of output columns. If set, only the `n_out_cols` first elements of the array column are output.
                             If <1, the number of output columns is equal to the length of the shortest array in `array_col`
    :param out_cols_prefix: Prefix for the output columns (uses `array_col` as the prefix unless set)
    :param bool out_1based_indexing: If set, the output column names indexes start at 1. Otherwise they start at 0.
    :return: dataframe with expanded columns
    :rtype: DataFrame
    """

    if out_cols_prefix is None:
        out_cols_prefix = array_col

    if num_out_cols < 1:
        num_out_cols = min([len(x) for x in df[array_col].values.tolist()])

    cols = ['{}{}'.format(out_cols_prefix, i + out_1based_indexing) for i in range(num_out_cols)]
    df[cols] = pd.DataFrame(df[array_col].values.tolist())[list(range(num_out_cols))]

    return df


def assign_population_pcs(
        pop_pca_scores: Union[hl.Table, pd.DataFrame],
        pc_cols: Union[hl.expr.ArrayExpression, List[str]],
        known_col: str = 'known_pop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[Union[hl.Table, pd.DataFrame], RandomForestClassifier]:
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.
    As input, this function can either take:
    - A Hail Table (typically the output of `hwe_normalized_pca`). In this case,
        - `pc_cols` should be an ArrayExpression of Floats where each element is one of the PCs to use.
        - A Hail Table will be returned as output
    - A Pandas DataFrame. In this case:
        - Each PC should be in a separate column and `pc_cols` is the list of all the columns containing the PCs to use.
        - A pandas DataFrame is returned as output
    Note
    ----
    If you have a Pandas Dataframe and have all PCs as an array in a single column, the
    `expand_pd_array_col` can be used to expand this column into multiple `PC` columns.
    :param Table or DataFrame pop_pc_pd: Input Hail Table or Pandas Dataframe
    :param ArrayExpression or list of str pc_cols: Columns storing the PCs to use
    :param str known_col: Column storing the known population labels
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Hail Table or Pandas Dataframe (depending on input) containing sample IDs and imputed population labels, trained random forest model
    :rtype: (Table or DataFrame, RandomForestClassifier)
    """

    hail_input = isinstance(pop_pca_scores, hl.Table)
    if hail_input:
        pop_pc_pd = pop_pca_scores.select(
            known_col,
            pca_scores=pc_cols
        ).to_pandas()
        pop_pc_pd = expand_pd_array_col(pop_pc_pd, 'pca_scores', out_cols_prefix='PC')
        pc_cols = [col for col in pop_pc_pd if col.startswith('PC')]
    else:
        pop_pc_pd = pop_pca_scores

    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    N = len(train_data)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].values
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    pop_pc_pd = pd.concat([pop_pc_pd, probs], axis=1)
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    pop_pc_pd = pop_pc_pd.drop(pc_cols, axis='columns')

    print("Found the following sample count after population assignment: {}".format(
        ", ".join(f'{pop}: {count}' for pop, count in Counter(pop_pc_pd[output_col]).items())
    ))

    if hail_input:
        pops_ht = hl.Table.from_pandas(pop_pc_pd, key=list(pop_pca_scores.key))
        pops_ht.annotate_globals(
            assign_pops_from_pc_params=hl.struct(
                min_assignment_prob=min_prob
            )
        )
        return pops_ht, pop_clf
    else:
        return pop_pc_pd, pop_clf


