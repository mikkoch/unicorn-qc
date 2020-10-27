import sys
import json
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


######################################################################################################
# Defining the list of input bfiles, output root directory                                           #
######################################################################################################
bfiles = ['gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_bel1/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final',
          'gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_bel2/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final',
          'gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_chop/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final',
          'gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_germ/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final',
          'gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_nidk/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final',
          'gs://pgc-crohns-data/qced_genotyping_data/cd_data_cd_wtcc/hg19/hg19_ready4QC_AF_HRC_sa/hg19.all.final']

output_root_directories = ["gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/bel1",
                           "gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/bel2",
                           "gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/chop",
                           "gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/germ",
                           "gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/nidk",
                           "gs://unicorn-qc/pre-imputation-qc/cohort-qc-crohns/wtcc"]
mt_1kg_path = "gs://unicorn-resources/1000_genomes/pop_4pop_mix_SEQ.mt"
mt_1kg_eur_path = "gs://unicorn-resources/1000_genomes/pop_euro_eur_with_aj.mt"



######################################################################################################
# Defining functions for cohort QC                                                                   #
######################################################################################################
def __main__(bfile: str, mt_1kg_path: str, mt_1kg_eur_path: str, output_root_directory: str,
             r2: float, pirsh: float, scrsh: float, scrsh_chr: float, vcrsh: float,
             mafrsh: float, hwersh: float, srh: int=50, n_partitions: int=200):
    """
        This function performs cohort-level QC including the following QC filters:
        - sample call rate > 0.98
        - sample call rate for each chromosome > 0.50
        - inbreeding coefficient should be less than 3 std from its mean
        - remove related samples (pi-hcat threshold = 0.0625)
        - remove sex-error samples
        - variant call rate > 0.98
        - minor allele frequencies > 0.01
        - within each population, p-hwe > 1e-4
        - autosome only
        This function uses a random forest classifier trained from 1000 genomes (AMR, AFR, EUR, EAS) to predict the continental ancestry for each sample.
        The detailed procedure includes:
            - filter to variants that are common in both dataset & LD prune
            - PCA with 1000 genomes (the input 1000 genomes is Stephan's cleaned version of 1KG
            - train random forest classifier with first 6 PCs
            - project genotyping data onto 1KG
            - predict the continental ancestry label using the projected PCs
        Within European population, uses another random forest classifier trained from 1000 genomes to predict whether they're from mainland European,
        Finland or Ashkenazi Jewish samples
        Note
        ----
        :param str bfile: path to input PLINK bfile
        :param str mt_1kg_path: path to 1000 genomes reference
        :param str mt_1kg_eur_path: path to 1000 genomes European reference
        :param str output_root_directory
        :param float r2: LD pruning threshold
        :param float pirsh: pihat threshold for filtering related samples
        :param float scrsh: sample call rate threshold
        :param float scrsh_chr: per chromosome sample call rate
        :param float vcrsh: variant call rate threshold
        :param float mafrsh: minor allele frequencies threshold
        :param float hwersh: hardy weinberg equilibrium test pvalues threshold
        :param int srh: cut-off for minimum sample size
        """

    mt = hl.import_plink(bed=bfile+'.bed', bim=bfile+'.bim', fam=bfile+'.fam', min_partitions=n_partitions)
    mt_1kg = hl.read_matrix_table(mt_1kg_path)
    mt_1kg_eur = hl.read_matrix_table(mt_1kg_eur_path)
    directory_structure = build_cohort_qc_directory_structure(output_root_directory)
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # If sample size is less than a threshold, skip
    if n_samples_before_qc < srh: return

    # Perform sample QC
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)
    mt = calculate_per_chr_scrt(mt=mt, scrt_col='scr_chr')  # calculate per-chromosome call rate for each sample
    mt = calculate_inbreeding_coefficients(mt=mt, ib_col='ib', ib_global='ib_stats')  # calculate inbreeding coefficient
    mt = mt.cache()

    # Calculate number of samples & variants that fail each filter
    n_scrsh     = mt.aggregate_cols(hl.agg.count_where(mt.sample_qc.call_rate < scrsh))
    n_scrsh_chr = mt.aggregate_cols(hl.agg.count_where(mt.scr_chr < scrsh_chr))
    n_ib        = mt.aggregate_cols(hl.agg.count_where(hl.abs(mt.ib.f_stat-mt.ib_stats.mean) > 3 * mt.ib_stats.stdev))
    n_vcrsh     = mt.aggregate_rows(hl.agg.count_where(mt.variant_qc.call_rate < vcrsh))
    n_mafsh     = mt.aggregate_rows(hl.agg.count_where((mt.variant_qc.AF[0] < mafrsh) | (mt.variant_qc.AF[1] < mafrsh)))

    # Histogram of QC metrics before QC
    histogram(mt.cols(),
              location='sample_qc.call_rate',
              plot_path=directory_structure["summary"]+"/hist_before_QC_sample_call_rate.png")
    histogram(mt.cols(),
              location='ib.f_stat',
              plot_path=directory_structure["summary"]+"/hist_before_QC_inbreeding_coefficient.png")
    histogram(mt.rows(),
              location='variant_qc.call_rate',
              plot_path=directory_structure["summary"]+"/hist_before_QC_variant_call_rate.png")

    # Filter samples & variants
    mt = mt.filter_cols((mt.sample_qc.call_rate > scrsh) &
                        (hl.abs(mt.ib.f_stat - mt.ib_stats.mean) <= 3 * mt.ib_stats.stdev) &
                        (mt.scr_chr > scrsh_chr), keep=True)
    mt = mt.filter_rows((mt.variant_qc.call_rate > vcrsh) &
                        (mt.variant_qc.AF[0] > mafrsh) &
                        (mt.variant_qc.AF[1] > mafrsh) &
                        (mt.locus.in_autosome()), keep=True)

    # Histogram of QC metrics after QC
    histogram(mt.cols(),
              location='sample_qc.call_rate',
              plot_path=directory_structure["summary"] + "/hist_after_QC_sample_call_rate.png")
    histogram(mt.cols(),
              location='ib.f_stat',
              plot_path=directory_structure["summary"] + "/hist_after_QC_inbreeding_coefficient.png")
    histogram(mt.rows(),
              location='variant_qc.call_rate',
              plot_path=directory_structure["summary"] + "/hist_after_QC_variant_call_rate.png")

    # LD pruning & Identify related samples
    pruned_variants_list = ld_prune(mt, r2=r2, pruned_variants_list=True)
    pruned_variants_list.write(directory_structure["pca"] + "/ld_pruned_variants.kt", overwrite=True)
    pruned_variants_list = hl.read_table(directory_structure["pca"] + "/ld_pruned_variants.kt")
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variants_list[mt.row_key]), keep=True)
    samples_to_remove = identify_related_samples(pruned_mt, pihat_threshold=pirsh)

    # Filter related samples
    n_related = samples_to_remove.count()
    mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]), keep=False)

    # Project onto 1000 genomes
    pruned_1kg_mt = mt_1kg.filter_rows(hl.is_defined(pruned_variants_list[mt_1kg.row_key]), keep=True)
    scores_1kg_ht, loadings_1kg_ht = pca(mt=pruned_1kg_mt, n_evecs=6, remove_outliers=False)
    scores_ht = pca_project(mt=mt, loadings_ht=loadings_1kg_ht, correct_shrinkage=True)
    scores_ht = scores_1kg_ht.union(scores_ht)
    scores_ht = scores_ht.annotate(super_population=mt_1kg.index_cols(scores_ht.key).super_population)

    # Use Random forest classifier to assign ancestry
    pops_ht, pop_clf = assign_population_pcs(
        pop_pca_scores=scores_ht,
        pc_cols=scores_ht.scores,
        known_col='super_population',
        min_prob=0.9,
        prop_train=0.8)

    # Plot ancestry assignment result
    pops_ht = pops_ht.annotate(PC1=pops_ht.pca_scores[0],
                               PC2=pops_ht.pca_scores[1],
                               PC3=pops_ht.pca_scores[2],
                               PC4=pops_ht.pca_scores[3],
                               PC5=pops_ht.pca_scores[4],
                               PC6=pops_ht.pca_scores[5])
    scatter(ht=pops_ht, x_location='PC1', y_location='PC2', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC1_PC2.png')
    scatter(ht=pops_ht, x_location='PC1', y_location='PC3', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC1_PC3.png')
    scatter(ht=pops_ht, x_location='PC2', y_location='PC3', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC2_PC3.png')
    scatter(ht=pops_ht, x_location='PC3', y_location='PC4', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC3_PC4.png')
    scatter(ht=pops_ht, x_location='PC4', y_location='PC5', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC4_PC5.png')
    scatter(ht=pops_ht, x_location='PC5', y_location='PC6', color_location='pop',
            plot_path=directory_structure['pca'] + '/SCATTER_1KG_PC5_PC6.png')
    pops_ht.write(directory_structure['pca'] + '/scores_1kg.kt', overwrite=True)
    mt = mt.annotate_cols(pop=pops_ht[mt.col_key].pop)

    # For samples assigning to EUR population, trying to project them onto 1KG EUR,
    # and classify into eur(mainland), FIN and AJ
    if mt.aggregate_cols(hl.agg.count_where(mt.pop == 'eur')) > 0:
        mt_eur = mt.filter_cols(mt.pop == 'eur', keep=True)
        mt_1kg_eur_mt = mt_1kg_eur.filter_rows(hl.is_defined(pruned_variants_list[mt_1kg_eur.row_key]), keep=True)
        scores_1kg_eur_ht, loadings_1kg_eur_ht = pca(mt=mt_1kg_eur_mt, n_evecs=6, remove_outliers=False)
        scores_ht = pca_project(mt=mt_eur, loadings_ht=loadings_1kg_eur_ht, correct_shrinkage=True)
        scores_ht = scores_1kg_eur_ht.union(scores_ht)
        scores_ht = scores_ht.annotate(population=mt_1kg_eur.index_cols(scores_ht.key).population)

        # Use Random forest classifier to assign ancestry
        pops_ht, pop_clf = assign_population_pcs(
            pop_pca_scores=scores_ht,
            pc_cols=scores_ht.scores,
            known_col='population',
            min_prob=0.9,
            prop_train=0.8)

        # Plot ancestry assignment result
        pops_ht = pops_ht.annotate(PC1=pops_ht.pca_scores[0],
                                   PC2=pops_ht.pca_scores[1],
                                   PC3=pops_ht.pca_scores[2],
                                   PC4=pops_ht.pca_scores[3],
                                   PC5=pops_ht.pca_scores[4],
                                   PC6=pops_ht.pca_scores[5])
        scatter(ht=pops_ht, x_location='PC1', y_location='PC2', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC1_PC2.png')
        scatter(ht=pops_ht, x_location='PC1', y_location='PC3', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC1_PC3.png')
        scatter(ht=pops_ht, x_location='PC2', y_location='PC3', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC2_PC3.png')
        scatter(ht=pops_ht, x_location='PC3', y_location='PC4', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC3_PC4.png')
        scatter(ht=pops_ht, x_location='PC4', y_location='PC5', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC4_PC5.png')
        scatter(ht=pops_ht, x_location='PC5', y_location='PC6', color_location='pop',
                plot_path=directory_structure['pca'] + '/SCATTER_1KG_EUR_PC5_PC6.png')
        pops_ht.write(directory_structure['pca']+'/scores_1kg_eur.kt', overwrite=True)
        mt = mt.transmute_cols(pop=hl.cond(mt.pop == 'eur', pops_ht[mt.col_key].pop, mt.pop))

    # Filter out samples that are not assigned to any populations
    mt = mt.filter_cols(mt.pop == 'oth', keep=False)

    # Count number of samples in each population
    pops = mt.aggregate_cols(hl.agg.counter(mt.pop))

    # Calculate p-values of Hardy Weinberg Equilibrium within each population
    for pop, count in pops.items():
        if count < srh:  # if number of samples in this population is lower than certain threshold, skip HWE test
            continue
        mt = mt.annotate_rows(**{"pop_"+pop:hl.agg.filter(mt.pop == pop, hl.agg.hardy_weinberg_test(mt.GT))})

    # Calculate minimum hwe p-values across populations
    mt = mt.annotate_rows(het_freq_hwe=hl.min([mt['pop_'+pop].het_freq_hwe for pop, count in pops.items() if count > srh]),
                               p_value=hl.min([mt['pop_'+pop].p_value for pop, count in pops.items() if count > srh]))

    # Count number of variants failing HWE filter
    n_hrsh = mt.aggregate_rows(hl.agg.count_where((mt.het_freq_hwe < hwersh) | (mt.p_value < hwersh)))

    # Filter out variants that fails HWE
    mt = mt.filter_rows((mt.het_freq_hwe > hwersh) & (mt.p_value > hwersh), keep=True)

    # Count samples & variants after QC
    n_variants_after_qc, n_samples_after_qc = mt.count()

    # create a dictionary storing the meta information
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc,
            "Sample QC Summary": {
                "IDs: call rate < %s" % scrsh: n_scrsh,
                "IDs: minimum per-chromosome call rate < %s" %scrsh_chr: n_scrsh_chr,
                "IDs: pi-hat > %s" % pirsh: n_related,
                "IDs: inbreeding coefficient is 3 stdev from its mean": n_ib},
            "Variant QC Summary": {
                "SNPs: call rate < %s" % vcrsh: n_vcrsh,
                "SNPs: minor allele frequency < %s" % mafrsh: n_mafsh,
                "SNPs: HWE p-values < %s" % hwersh: n_hrsh},
            "Population Assignment": {
                "EUR (mainland)": pops.get("eur(mainland)", 0),
                "Ashkenazi Jew": pops.get("aj", 0),
                "FIN": pops.get("fin", 0),
                "EAS": pops.get("asn", 0),
                "AFR": pops.get("afr", 0),
                "AMR": pops.get("amr", 0)
            }}
    mt.write(directory_structure["share"]+"/all.mt", overwrite=True)

    # Export meta
    with hl.hadoop_open(directory_structure["summary"] + "/meta.json", 'w') as outfile:
        json.dump(meta, outfile)

    # Write hail MatrixTable
    for pop, count in pops.items():
        pop_mt = mt.filter_cols(mt.pop == pop)
        pop_mt.write(directory_structure["share"]+"/{pop}.mt".format(pop=pop), overwrite=True)




######################################################################################################
# Run cohort QC pipeline                                                                             #
######################################################################################################
input_pandas = pd.DataFrame({'bfile':bfiles, 'out_directory': output_root_directories})
for index, row in input_pandas.iterrows():
    __main__(bfile=row['bfile'],
             mt_1kg_path=mt_1kg_path,
             mt_1kg_eur_path=mt_1kg_eur_path,
             output_root_directory=row['out_directory'],
             r2=r2,
             pirsh=pihat_threshold,
             scrsh=sample_call_rate_threshold,
             scrsh_chr=call_rate_per_chromosome_threshold,
             vcrsh=variant_call_rate_threshold,
             mafrsh=maf_threshold,
             hwersh=hwe_threshold)