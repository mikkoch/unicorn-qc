import hail as hl

def groupby_populations(
    mt: hl.MatrixTable,
    cluster_ht: hl.Table,
    rows_cluster_name:str='cluster',
    output_rows_cluster_name: str= 'cluster') ->hl.MatrixTable:

    mt = mt.annotate_cols(**{
        output_rows_cluster_name: cluster_ht[mt.col_key][rows_cluster_name]})
    mt = mt.filter_cols(
        hl.is_defined(mt[output_rows_cluster_name]))
    pops_mt = (mt
        .group_cols_by(mt[output_rows_cluster_name])
        .aggregate(
            n_alt_alleles=hl.agg.sum(mt.GT.n_alt_alleles()),
            n_called=hl.agg.count_where(hl.is_defined(mt.GT))))

    pops_mt = pops_mt.annotate_entries(
        n_ref_alleles=2*pops_mt.n_called-pops_mt.n_alt_alleles)

    pops_mt = pops_mt.annotate_entries(
        AF=pops_mt.n_alt_alleles / (2.0*pops_mt.n_called))

    pops_mt = pops_mt.annotate_entries(
        AF_VAR=pops_mt.AF*(1-pops_mt.AF)/ (2*pops_mt.n_called),
        AF_SQR=pops_mt.AF**2)

    return pops_mt


def cmh_test(
    mafd_co:hl.MatrixTable,
    mafd_ca:hl.MatrixTable) -> hl.Table:

    mafd_co = mafd_co.filter_cols(
        hl.is_defined(mafd_ca.index_cols(mafd_co.col_key)),
        keep=True)
    mafd_co = mafd_co.filter_rows(
        hl.is_defined(mafd_ca.index_rows(mafd_co.row_key)),
        keep=True)

    mafd = mafd_co.annotate_entries(
        n_alt_gt_case=mafd_ca[mafd_co.row_key, mafd_co.col_key].n_alt_alleles,
        n_alt_gt_cont=mafd_co.n_alt_alleles,
        n_ref_gt_case=mafd_ca[mafd_co.row_key, mafd_co.col_key].n_ref_alleles,
        n_ref_gt_cont=mafd_co.n_ref_alleles)

    mafd = mafd.annotate_entries(
        m_1i=mafd.n_alt_gt_case+mafd.n_alt_gt_cont,
        m_2i=mafd.n_ref_gt_case+mafd.n_ref_gt_cont,
        n_1i=mafd.n_alt_gt_case+mafd.n_ref_gt_case,
        n_2i=mafd.n_alt_gt_cont+mafd.n_ref_gt_cont,
        t_i =mafd.n_alt_gt_case+mafd.n_alt_gt_cont+mafd.n_ref_gt_case+mafd.n_ref_gt_cont)

    mafd = mafd.filter_entries(
        (mafd.n_alt_gt_case==0) &
        (mafd.n_alt_gt_cont==0), keep=False)

    mafd = mafd.annotate_rows(
        numer=hl.agg.sum(mafd.n_alt_gt_case - mafd.n_1i*mafd.m_1i/mafd.t_i),
        denom=hl.agg.sum(mafd.n_1i*mafd.n_2i*mafd.m_1i*mafd.m_2i/(mafd.t_i*mafd.t_i*(mafd.t_i-1))))

    mafd = mafd.annotate_rows(chi2_stat=mafd.numer**2 / mafd.denom)

    chi2_stats_ht = mafd.rows().select('chi2_stat')

    return chi2_stats_ht



def inference_frequentist(
        ca_mt: hl.MatrixTable,
        pops_mt: hl.MatrixTable,
        probs_ht: hl.Table,
        rows_cluster_name:str='cluster'):

    populations = pops_mt[rows_cluster_name].collect()
    probs_mt = probs_ht.to_matrix_table_row_major(populations, entry_field_name='prob')

    probs_mt = probs_mt.annotate_cols(
        prob_sum=hl.agg.sum(probs_mt.prob),
        prob_sqr_sum=hl.agg.sum(probs_mt.prob ** 2),
        prob_sum_sqr=hl.agg.sum(probs_mt.prob)** 2)

    pops_mt = pops_mt.annotate_cols(
        prob_sum=probs_mt.index_cols(pops_mt.col_key).prob_sum,
        prob_sqr_sum=probs_mt.index_cols(pops_mt.col_key).prob_sqr_sum,
        prob_sum_sqr=probs_mt.index_cols(pops_mt.col_key).prob_sum_sqr)


    pops_mt = pops_mt.annotate_rows(
        y_sum_exp=2 * hl.agg.sum(pops_mt.prob_sum * pops_mt.AF),
        y_sum_var=hl.agg.sum(
            2 * pops_mt.AF * pops_mt.prob_sum -
            2 * pops_mt.AF_SQR * pops_mt.prob_sqr_sum +
            4 * pops_mt.AF_VAR * pops_mt.prob_sum_sqr))

    ca_mt = ca_mt.annotate_rows(
        y_sum=hl.agg.sum(ca_mt.GT.n_alt_alleles()),
        y_sum_exp=pops_mt.index_rows(ca_mt.row_key).y_sum_exp,
        y_sum_var=pops_mt.index_rows(ca_mt.row_key).y_sum_var)

    ca_mt = ca_mt.annotate_rows(
        chi2_stat=((ca_mt.y_sum - ca_mt.y_sum_exp)**2) / ca_mt.y_sum_var)
    chi2_stats_ht = ca_mt.rows().select('chi2_stat')

    return chi2_stats_ht