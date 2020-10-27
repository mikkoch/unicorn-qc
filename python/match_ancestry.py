import sys
import json

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *


def match_ancestry(path_to_mt: str, mt_1kg_path: str, mt_1kg_eur_path: str, output_root_directory: str,
             r2: float, pihat: float, hwe_thresh: float, n_partitions: int = 200):
    """  
        This function uses a random forest classifier trained from 1000 genomes to predict the continental ancestry for each sample.
        Within European population, uses another random forest classifier trained from 1000 genomes to predict whether they're from mainland European,
        Finland or Ashkenazi Jewish samples
        Note
        ----
        :param str path_to_mt: path to matrixtable
        :param str mt_1kg_path: path to 1000 genomes reference
        :param str mt_1kg_eur_path: path to 1000 genomes European reference
        :param str output_root_directory
        :param float r2: LD pruning threshold
        :param float hwe_thresh: hardy weinberg equilibrium test pvalues threshold
        :param n_partitions: number of partitions
        """

    mt = hl.read_matrix_table(path_to_mt)
    mt_1kg = hl.read_matrix_table(mt_1kg_path)
    mt_1kg_eur = hl.read_matrix_table(mt_1kg_eur_path)
    directory_structure = build_cohort_qc_directory_structure(output_root_directory)
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # LD pruning & Identify related samples
    pruned_variants_list = ld_prune(mt, r2=r2, pruned_variants_list=True)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_variants_list[mt.row_key]), keep=True)
    samples_to_remove = identify_related_samples(pruned_mt, pihat_threshold=pihat)

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
        pop_pca_scores = scores_ht, pc_cols = scores_ht.scores,
        known_col = 'super_population', min_prob = 0.9, prop_train = 0.8)

    # Plot ancestry assignment result
    pops_ht = pops_ht.annotate(PC1=pops_ht.pca_scores[0],
                               PC2=pops_ht.pca_scores[1],
                               PC3=pops_ht.pca_scores[2],
                               PC4=pops_ht.pca_scores[3],
                               PC5=pops_ht.pca_scores[4],
                               PC6=pops_ht.pca_scores[5])
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
            pop_pca_scores=scores_ht, pc_cols = scores_ht.scores,
            known_col='population', min_prob = 0.9, prop_train = 0.8)

        pops_ht.write(directory_structure['pca']+'/scores_1kg_eur.kt', overwrite=True)
        mt = mt.transmute_cols(pop=hl.cond(mt.pop == 'eur', pops_ht[mt.col_key].pop, mt.pop))

    # Filter out samples that are not assigned to any populations
    mt = mt.filter_cols(mt.pop == 'oth', keep=False)

    # Count number of samples in each population
    pops = mt.aggregate_cols(hl.agg.counter(mt.pop))

    # Calculate p-values of Hardy Weinberg Equilibrium within each population
    for pop, count in pops.items():
        mt = mt.annotate_rows(**{"pop_"+pop:hl.agg.filter(mt.pop == pop, hl.agg.hardy_weinberg_test(mt.GT))})

    # Calculate minimum hwe p-values across populations
    mt = mt.annotate_rows(het_freq_hwe=hl.min([mt['pop_'+pop].het_freq_hwe for pop, count in pops.items()]),
                          p_value=hl.min([mt['pop_'+pop].p_value for pop, count in pops.items()]))

    # Count number of variants failing HWE filter
    n_hwe_variants = mt.aggregate_rows(hl.agg.count_where((mt.het_freq_hwe < hwe_thresh) | (mt.p_value < hwe_thresh)))

    # Filter out variants that fails HWE
    mt = mt.filter_rows((mt.het_freq_hwe > hwe_thresh) & (mt.p_value > hwe_thresh), keep = True)

    # Count samples & variants after QC
    n_variants_after_qc, n_samples_after_qc = mt.count()

    # create a dictionary storing the meta information
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc,
            "SNPs: HWE p-values < %s" % hwe_thresh: n_hwe_variants,
            "Population Assignment": {
                "EUR (mainland)": pops.get("eur(mainland)", 0),
                "Ashkenazi Jew": pops.get("aj", 0),
                "FIN": pops.get("fin", 0),
                "EAS": pops.get("asn", 0),
                "AFR": pops.get("afr", 0),
                "AMR": pops.get("amr", 0),
                "SAS": pops.get("sas", 0)
            }}
    mt.write(directory_structure["share"]+"/all.mt", overwrite = True)

    # Export meta
    with hl.hadoop_open(directory_structure["summary"] + "/match_ancestry_summary.json", 'w') as outfile:
        json.dump(meta, outfile)

    # Write hail MatrixTable
    for pop, count in pops.items():
        pop_mt = mt.filter_cols(mt.pop == pop)
        pop_mt.write(directory_structure["share"]+"/{pop}.mt".format(pop=pop), overwrite=True)


