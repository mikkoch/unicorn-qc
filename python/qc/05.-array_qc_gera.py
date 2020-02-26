import sys
import json
from typing import *
from functools import reduce
sys.path.insert(0, '/home/danfengc/unicorn')

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *


######################################################################################################
# Defining QC parameters                                                                             #
######################################################################################################
r2 = 0.2  # LD pruning threshold
pihat_threshold = 0.0625  # related/duplicated samples threshold (plink --genome)
sample_call_rate_threshold = 0.98  # sample call rate
call_rate_per_chromosome_threshold = 0.50  # per chromosome sample call rate (Michigan Imputation Server)
variant_call_rate_threshold = 0.98  # variant call rate
maf_threshold = 0.01  # MAF threshold
hwe_threshold = 1e-4  # Hardy Weinberg threshold
aso_threshold = 1e-4  # 1KG association threshold


mt_1kg_eur_path = "gs://unicorn-resources/1000_genomes/pop_euro_eur_SEQ.mt"
output_root_directory = "gs://unicorn-qc/pre-imputation-qc/array-qc"



def __main__(mt_path: str, cohort_label: str, output_root_directory:str, population:str, array:str,
             mt_1kg_eur_path:str, r2:float, pirsh: float, scrsh:float, scrsh_chr: float, vcrsh: float, mafrsh: float,
             hwersh: float, assrsh: float):
    """
        This function performs cohort-level QC including the following QC filters:
        - sample call rate > 0.98
        - sample call rate for each chromosome > 0.50
        - remove related samples (pi-hat threshold = 0.0625)
        - variant call rate > 0.98
        - minor allele frequencies > 0.01
        - within each population, p-hwe > 1e-4
        - pseudo case-control GWAS (tagging samples from one cohort as cases, from the other cohorts as controls) p-values > 1e-4
        - pseudo case-control GWAS against 1KG p-values > 1e-4 (currently it is EUR only)
        Note
        In order to convert VCF to zipped VCF, run the following commands:
            sudo apt-get update
            sudo apt-get install vcftools -y
            sudo apt-get install tabix -y
        ----
        :param List[str] mt_path_list: a list of path to MatrixTable that are from same population & same array
        :param List[str] cohort_labels: a list of cohort labels
        :param str output_root_directory
        :param str population
        :param str array
        :param str mt_1kg_eur_path: path to 1KG EUR population
        :param float r2: LD pruning threshold
        :param float pirsh: pihat threshold for filtering related samples
        :param float scrsh: sample call rate threshold
        :param float scrsh_chr: per chromosome sample call rate
        :param float vcrsh: variant call rate threshold
        :param float mafrsh: minor allele frequencies threshold
        :param float hwersh: hardy weinberg equilibrium test pvalues threshold
        :param float assrsh: p-value threshold for pseudo GWAS
        """

    # Merge MatrixTable from same population and same genotyping array
    mt = hl.read_matrix_table(mt_path)
    mt = mt.annotate_cols(cohort=cohort_label)

    # Read 1KG EUR mainland matrix table
    mt_1kg_eur = hl.read_matrix_table(mt_1kg_eur_path)
    mt_1kg_eur = mt_1kg_eur.filter_cols(mt_1kg_eur.population == 'eur(mainland)')

    # Build Directory structure
    directory_structure = build_array_qc_directory_structure(output_root_directory, population=population, array=array)

    # Count number of samples and variants before QC
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # Perform sample QC
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)
    mt = calculate_per_chr_scrt(mt=mt, scrt_col='scr_chr')

    # Calculate number of samples & variants that fail each filter
    n_scrsh     = mt.aggregate_cols(hl.agg.count_where(mt.sample_qc.call_rate < scrsh))
    n_scrsh_chr = mt.aggregate_cols(hl.agg.count_where(mt.scr_chr < scrsh_chr))
    n_vcrsh     = mt.aggregate_rows(hl.agg.count_where(mt.variant_qc.call_rate < vcrsh))
    n_mafsh     = mt.aggregate_rows(hl.agg.count_where((mt.variant_qc.AF[0] < mafrsh) | (mt.variant_qc.AF[1] < mafrsh)))
    n_hwesh     = mt.aggregate_rows(hl.agg.count_where((mt.het_freq_hwe < hwersh) | (mt.p_value < hwersh)))

    # Filter samples & variants
    mt = mt.filter_cols(
        (mt.sample_qc.call_rate > scrsh) &
        (mt.scr_chr > scrsh_chr), keep=True)
    mt = mt.filter_rows(
        (mt.variant_qc.call_rate > vcrsh) &
        (mt.variant_qc.AF[0] > mafrsh) &
        (mt.variant_qc.AF[1] > mafrsh) &
        (mt.het_freq_hwe > hwersh) &
        (mt.p_value > hwersh), keep=True)

    # LD pruning
    pruned_variants_list = hl.read_table("gs://unicorn-qc/pre-imputation-qc/cohort-qc/gera_eur_10000/pca/ld_pruned_variants.kt")
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variants_list[mt.row_key]), keep=True)
    #samples_to_remove = identify_related_samples(pruned_mt, pihat_threshold=pirsh)

    # Filter related samples
    n_related = 0

    # PCA
    scores_ht, _ = pca(mt=pruned_mt, n_evecs=6, remove_outliers=False)
    scores_ht.write(directory_structure["pca"] + "/scores.ht", overwrite=True)

    # Plot PCs with cohort label
    scores_ht = scores_ht.annotate(cohort=mt.index_cols(scores_ht.key).cohort)
    scores_ht = scores_ht.annotate(
        PC1=scores_ht.scores[0],
        PC2=scores_ht.scores[1],
        PC3=scores_ht.scores[2],
        PC4=scores_ht.scores[3],
        PC5=scores_ht.scores[4],
        PC6=scores_ht.scores[5])
    scatter(ht=scores_ht, x_location='PC1', y_location='PC2', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC1_PC2_BEFORE_QC.png')
    scatter(ht=scores_ht, x_location='PC1', y_location='PC3', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC1_PC3_BEFORE_QC.png')
    scatter(ht=scores_ht, x_location='PC2', y_location='PC3', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC2_PC3_BEFORE_QC.png')
    scatter(ht=scores_ht, x_location='PC3', y_location='PC4', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC3_PC4_BEFORE_QC.png')
    scatter(ht=scores_ht, x_location='PC4', y_location='PC5', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC4_PC5_BEFORE_QC.png')
    scatter(ht=scores_ht, x_location='PC5', y_location='PC6', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC5_PC6_BEFORE_QC.png')


    # Pseudo Case-Control Analysis against 1KG EUR
    n_cs1kg = 0
    if population == "eur_mainland":
        pruned_1kg_eur_mt = mt_1kg_eur.filter_rows(hl.is_defined(pruned_variants_list[mt_1kg_eur.row_key]), keep=True)
        scores_1kg_eur_ht, loadings_1kg_eur_ht = pca(mt=pruned_1kg_eur_mt, n_evecs=6, remove_outliers=False)
        mt_merged = mt.select_cols().select_rows().union_cols(mt_1kg_eur.select_cols().select_rows())
        scores_ht = scores_ht.select('scores').union(scores_1kg_eur_ht)

        mt_merged = mt_merged.annotate_cols(scores=scores_ht[mt_merged.col_key].scores, is_case=hl.cond(hl.is_defined(mt.index_cols(mt_merged.col_key)), 1, 0))
        result_ht = hl.linear_regression_rows(y=mt_merged.is_case, x=mt_merged.GT.n_alt_alleles(), covariates=[1, mt_merged.scores[0], mt_merged.scores[1], mt_merged.scores[2], mt_merged.scores[3], mt_merged.scores[4], mt_merged.scores[5]])

        # Annotate p-values from comparison against 1KG EUR
        mt = mt.annotate_rows(pvalues_1kg_assoc=result_ht[mt.row_key].p_value)

        # Count number of variants failing comparison against 1KG EUR
        n_cs1kg = mt.aggregate_rows(hl.agg.count_where(mt.pvalues_1kg_assoc < assrsh))

        # Filter variants
        mt = mt.filter_rows((mt.pvalues_1kg_assoc > assrsh) |
                            (hl.is_missing(mt.pvalues_1kg_assoc)), keep=True)

        # PLOT projecting onto 1KG EUR
        scores_ht = scores_ht.annotate(PC1=scores_ht.scores[0],
                                       PC2=scores_ht.scores[1],
                                       PC3=scores_ht.scores[2],
                                       PC4=scores_ht.scores[3],
                                       PC5=scores_ht.scores[4],
                                       PC6=scores_ht.scores[5])
        scatter(ht=scores_ht, x_location='PC1', y_location='PC2',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC1_PC2.png')
        scatter(ht=scores_ht, x_location='PC1', y_location='PC3',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC1_PC3.png')
        scatter(ht=scores_ht, x_location='PC2', y_location='PC3',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC2_PC3.png')
        scatter(ht=scores_ht, x_location='PC3', y_location='PC4',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC3_PC4.png')
        scatter(ht=scores_ht, x_location='PC4', y_location='PC5',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC4_PC5.png')
        scatter(ht=scores_ht, x_location='PC5', y_location='PC6',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC5_PC6.png')
        scores_ht.write(directory_structure['pca'] + '/scores_1kg.kt', overwrite=True)

    # LD pruning
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variants_list[mt.row_key]), keep=True)

    # PCA
    scores_ht, _ = pca(mt=pruned_mt, n_evecs=6, remove_outliers=True, sigma_thresh=5, n_outlieriters=1)
    scores_ht.write(directory_structure["pca"] + "/scores.ht", overwrite=True)

    # Plot PCs with cohort label
    scores_ht = scores_ht.annotate(cohort=mt.index_cols(scores_ht.key).cohort)
    scores_ht = scores_ht.annotate(PC1=scores_ht.scores[0],
                                   PC2=scores_ht.scores[1],
                                   PC3=scores_ht.scores[2],
                                   PC4=scores_ht.scores[3],
                                   PC5=scores_ht.scores[4],
                                   PC6=scores_ht.scores[5])
    scatter(ht=scores_ht, x_location='PC1', y_location='PC2', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC1_PC2_AFTER_QC.png')
    scatter(ht=scores_ht, x_location='PC1', y_location='PC3', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC1_PC3_AFTER_QC.png')
    scatter(ht=scores_ht, x_location='PC2', y_location='PC3', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC2_PC3_AFTER_QC.png')
    scatter(ht=scores_ht, x_location='PC3', y_location='PC4', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC3_PC4_AFTER_QC.png')
    scatter(ht=scores_ht, x_location='PC4', y_location='PC5', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC4_PC5_AFTER_QC.png')
    scatter(ht=scores_ht, x_location='PC5', y_location='PC6', color_location='cohort',
            plot_path=directory_structure['pca'] + '/SCATTER_PC5_PC6_AFTER_QC.png')

    # Count PCA outliers
    n_pca_outliers = mt.aggregate_cols(hl.agg.count_where(~hl.is_defined(scores_ht[mt.col_key])))

    # Filter PCA outliers
    mt = mt.filter_cols(hl.is_defined(scores_ht[mt.col_key]), keep=True)

    # Count number of samples and variants after QC
    n_variants_after_qc, n_samples_after_qc = mt.count()

    # Write QC'ed data to google bucket
    mt.write(directory_structure["share"]+"/final.mt", overwrite=True)
    #mt = hl.read_matrix_table(directory_structure["share"]+"/final.mt")

    # Write QC meta information
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc,
            "Genotyping Array": array,
            "Population": population,
            "Sample QC Summary": {
                "IDs: call rate < %s" % scrsh: n_scrsh,
                "IDs: minimum per-chromosome call rate < %s" % scrsh_chr: n_scrsh_chr,
                "IDs: pi-hat > %s" % pirsh: n_related,
                "IDs: PCA outliers": n_pca_outliers},
            "Variant QC Summary": {
                "SNPs: call rate < %s" % vcrsh: n_vcrsh,
                "SNPs: minor allele frequency < %s" % mafrsh: n_mafsh,
                "SNPs: HWE p-values < %s" % hwersh: n_hwesh,
                "SNPs: Cross cohorts pseudo case-control across p-value < %s" % assrsh: n_cs1kg}}
    with hl.hadoop_open(directory_structure["summary"] + "/meta.json", 'w') as outfile:
        json.dump(meta, outfile)

    # Convert to VCF files
    convert_to_vcf(mt=mt,
                   output_root_for_unzip_vcf=directory_structure['vcf']+"/unzip",
                   output_root_for_zip_vcf=directory_structure['vcf']+"/zip",
                   basename='final')


######################################################################################################
# Run QC pipeline                                                                                    #
######################################################################################################

__main__(
    mt_path="gs://unicorn-qc/pre-imputation-qc/cohort-qc/gera_eur/share/eur(mainland).mt/",
    cohort_label='gera_eur',
    output_root_directory="gs://unicorn-qc/pre-imputation-qc/array-qc",
    population="eur_mainland",
    array="Axiom_KP_UCSF_EUR",
    mt_1kg_eur_path="gs://unicorn-resources/1000_genomes/pop_euro_eur_SEQ.mt",
    r2=r2,
    pirsh=pihat_threshold,
    scrsh=sample_call_rate_threshold,
    scrsh_chr=call_rate_per_chromosome_threshold,
    vcrsh=variant_call_rate_threshold,
    mafrsh=maf_threshold,
    hwersh=hwe_threshold,
    assrsh=aso_threshold)

