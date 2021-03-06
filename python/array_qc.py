import sys
import json
from typing import *
from functools import reduce

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *


def array_qc(path_to_mt: str, output_root_directory: str, mt_1kg_eur_path: str, r2: float, 
    pihat: float, sample_call_rate_thresh:float, sample_call_rate_thresh_chr: float, 
    variant_call_rate_thresh: float, maf_thresh: float, hwe_thresh: float, ass_thresh: float):
    """
        This function performs array-level QC including the following QC filters:
        - sample call rate
        - sample call rate for each chromosome 
        - variant call rate 
        - minor allele frequencies 
        - hardy weinberg equilibrium
        - pseudo case-control GWAS (tagging samples from one cohort as cases, from the other cohorts as controls)
        - pseudo case-control GWAS against 1KG 
        ----
        :param path_to_str: path to martrixtable
        :param str output_root_directory
        :param str mt_1kg_eur_path: path to 1KG EUR population
        :param float r2: LD pruning threshold
        :param float pihat: pihat threshold for filtering related samples
        :param float sample_call_rate_thresh: sample call rate threshold
        :param float sample_call_rate_thresh_chr: per chromosome sample call rate
        :param float variant_call_rate_thresh: variant call rate threshold
        :param float maf_thresh: minor allele frequencies threshold
        :param float hwe_thresh: hardy weinberg equilibrium test pvalues threshold
        :param float ass_thresh: p-value threshold for pseudo GWAS
        """

    mt = hl.read_matrix_table(path_to_mt)
    mt_1kg_eur = hl.read_matrix_table(mt_1kg_eur_path)
    mt_1kg_eur = mt_1kg_eur.filter_cols(mt_1kg_eur.population == 'eur(mainland)')

    # Build Directory structure
    directory_structure = build_array_qc_directory_structure(output_root_directory, population=population, array=array)

    # Count number of samples and variants before QC
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # Perform sample QC
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)
    mt = calculate_per_chr_sample_call_rate(mt = mt, sample_call_rate_col = 'scr_chr') 
    mt = mt.cache()

    # Calculate number of samples & variants that fail each filter
    n_sample_call_rate_thresh     = mt.aggregate_cols(hl.agg.count_where(mt.sample_qc.call_rate < sample_call_rate_thresh))
    n_sample_call_rate_thresh_chr = mt.aggregate_cols(hl.agg.count_where(mt.scr_chr < sample_call_rate_thresh_chr))
    n_variant_call_rate_thresh     = mt.aggregate_rows(hl.agg.count_where(mt.variant_qc.call_rate < variant_call_rate_thresh))
    n_mafsh     = mt.aggregate_rows(hl.agg.count_where((mt.variant_qc.AF[0] < maf_thresh) | (mt.variant_qc.AF[1] < maf_thresh)))
    n_hwesh     = mt.aggregate_rows(hl.agg.count_where((mt.het_freq_hwe < hwe_thresh) | (mt.p_value < hwe_thresh)))

    # Filter samples & variants
    mt = mt.filter_cols((mt.sample_qc.call_rate > sample_call_rate_thresh) &
                        (mt.scr_chr > sample_call_rate_thresh_chr), keep=True)
    mt = mt.filter_rows((mt.variant_qc.call_rate > variant_call_rate_thresh) &
                        (mt.variant_qc.AF[0] > maf_thresh) &
                        (mt.variant_qc.AF[1] > maf_thresh) &
                        (mt.het_freq_hwe > hwe_thresh) &
                        (mt.p_value > hwe_thresh), keep=True)

    # LD pruning
    pruned_variants_list = ld_prune(mt, r2=r2, pruned_variants_list=True)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variants_list[mt.row_key]), keep=True)
    samples_to_remove = identify_related_samples(pruned_mt, pihat_threshold=pihat)

    # Filter related samples
    n_related = samples_to_remove.count()
    mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]), keep=False)

    # PCA
    scores_ht, _ = pca(mt=pruned_mt, n_evecs=6, remove_outliers=False)
    scores_ht.write(directory_structure["pca"] + "/scores.ht", overwrite=True)

    # Pseudo Case-Control Analysis
    n_pseudo_case_ontrol = 0
    pruned_1kg_eur_mt = mt_1kg_eur.filter_rows(hl.is_defined(pruned_variants_list[mt_1kg_eur.row_key]), keep=True)
    scores_1kg_eur_ht, loadings_1kg_eur_ht = pca(mt=pruned_1kg_eur_mt, n_evecs=6, remove_outliers=False)
    scores_ht = pca_project(mt=mt, loadings_ht=loadings_1kg_eur_ht, correct_shrinkage=True)
    mt = mt.annotate_cols(scores=scores_ht[mt.col_key].scores)
    if len(cohort_labels) > 1:
        for cohort in cohort_labels:
            mt = mt.annotate_cols(is_case=hl.cond(mt.cohort == cohort, 1, 0))
            result_ht = hl.linear_regression_rows(y=mt.is_case, x=mt.GT.n_alt_alleles(), covariates=[1, mt.scores[0], mt.scores[1], mt.scores[2], mt.scores[3], mt.scores[4], mt.scores[5]])
            mt = mt.annotate_rows(**{'pvalue_' + cohort: result_ht[mt.row_key].p_value})
            mt = mt.drop('is_case')

        # Calculate the minimum pvalues cross each each cross cohort comparison
        mt = mt.annotate_rows(pvalues_assoc=hl.min([mt['pvalue_' + str(cohort)] for cohort in cohort_labels]))

        # Count number of variants failing cross-cohort comparison
        n_pseudo_case_ontrol = mt.aggregate_rows(hl.agg.count_where(mt.pvalues_assoc < ass_thresh))

        # Filter variants
        mt = mt.filter_rows(mt.pvalues_assoc > ass_thresh, keep=True)

    # Pseudo Case-Control Analysis against 1KG EUR
    n_pseudo_1kg_case_control = 0
    if population == "eur_mainland":
        mt_merged = mt.select_cols().select_rows().union_cols(mt_1kg_eur.select_cols().select_rows())
        scores_ht = scores_ht.union(scores_1kg_eur_ht)

        mt_merged = mt_merged.annotate_cols(scores=scores_ht[mt_merged.col_key].scores, is_case=hl.cond(hl.is_defined(mt.index_cols(mt_merged.col_key)), 1, 0))
        result_ht = hl.linear_regression_rows(y=mt_merged.is_case, x=mt_merged.GT.n_alt_alleles(), covariates=[1, mt_merged.scores[0], mt_merged.scores[1], mt_merged.scores[2], mt_merged.scores[3], mt_merged.scores[4], mt_merged.scores[5]])

        # Annotate p-values from comparison against 1KG EUR
        mt = mt.annotate_rows(pvalues_1kg_assoc=result_ht[mt.row_key].p_value)

        # Count number of variants failing comparison against 1KG EUR
        n_pseudo_1kg_case_control = mt.aggregate_rows(hl.agg.count_where(mt.pvalues_1kg_assoc < ass_thresh))

        # Filter variants
        mt = mt.filter_rows((mt.pvalues_1kg_assoc > ass_thresh) |
                            (hl.is_missing(mt.pvalues_1kg_assoc)), keep=True)


    # Count number of samples and variants after QC
    n_variants_after_qc, n_samples_after_qc = mt.count()

    # Write QC'ed data to google bucket
    mt.write(directory_structure["share"]+"/final.mt", overwrite=True)

    # Write QC meta information
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc,
            "Genotyping Array": array,
            "Population": population,
            "Sample QC Summary": {
                "IDs: call rate < %s" % sample_call_rate_thresh: n_sample_call_rate_thresh,
                "IDs: minimum per-chromosome call rate < %s" % sample_call_rate_thresh_chr: n_sample_call_rate_thresh_chr,
                "IDs: pi-hat > %s" % pihat: n_related}, 
            "Variant QC Summary": {
                "SNPs: call rate < %s" % variant_call_rate_thresh: n_variant_call_rate_thresh,
                "SNPs: minor allele frequency < %s" % maf_thresh: n_mafsh,
                "SNPs: HWE p-values < %s" % hwe_thresh: n_hwesh,
                "SNPs: Cross cohorts pseudo case-control across p-value < %s" % ass_thresh: n_pseudo_case_ontrol + n_pseudo_1kg_case_control}}
    with hl.hadoop_open(directory_structure["summary"] + "/meta.json", 'w') as outfile:
        json.dump(meta, outfile)

