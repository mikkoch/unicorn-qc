import os
import sys
from typing import *

sys.path.insert(0, '/home/danfengc/unicorn')
from main.pca import *
from qc.utils import *
from qc.plotting import *

def __main__(
        path_to_mt: str,
        path_to_pruned_mt: str,
        path_to_scores_ht: str,
        path_to_loadings_ht: str,
        path_to_plots: str,
        R2_thresh_for_pca: float = 0.98,
        af_thresh_for_pca: float = 0.05,
        n_partitions_for_pruning: int = 100):

    # Read MatrixTable from Google bucket
    mt = hl.read_matrix_table(path_to_mt)

    # arrays
    arrs = [
        'info_Axiom_KP_UCSF_EUR',
        'info_Affymetrix_6.0',
        'info_Human550',
        'info_Human610660',
        'info_HumanOmni']
 
    # loop through arrays
    for arr in arrs:
        mt = mt.filter_rows(
            (mt[arr].R2 > R2_thresh_for_pca) &
            (mt[arr].MAF > af_thresh_for_pca), keep=True)

    # Count number of variants Left
    print('Number of variants left for MAF threshold: {af_thresh}; ' \
          'R2 threshold: {r2_thresh}; '.format(
        af_thresh=af_thresh_for_pca,
        r2_thresh=R2_thresh_for_pca))

    # Re-partition
    mt = mt.repartition(n_partitions=n_partitions_for_pruning)

    # LD pruning
    pruned_mt = ld_prune(mt=mt)

    # Export pruned MatrixTable to Google BUCKET
    pruned_mt.write(path_to_pruned_mt, overwrite=True)

    # PCA
    eigens, scores_ht, loadings_ht = compute_full_spectrum(
        mt = pruned_mt,
        n_evecs = 20, 
        remove_outliers = True, 
        n_outlieriters = 1, 
        sigma_thresh = 5)

    # annotate
    scores_ht = scores_ht.annotate_globals(
        full_eigens = eigens)
    loadings_ht = loadings_ht.annotate_globals(
        full_eigens = eigens)

    # Export scores & loadings
    scores_ht.write(path_to_scores_ht, overwrite=True)
    loadings_ht.write(path_to_loadings_ht, overwrite=True)

    # PCA scatter-plot
    scores_ht = scores_ht.annotate(
        PC1=scores_ht.scores[0],
        PC2=scores_ht.scores[1],
        PC3=scores_ht.scores[2],
        PC4=scores_ht.scores[3],
        PC5=scores_ht.scores[4],
        PC6=scores_ht.scores[5],
        PC7=scores_ht.scores[6],
        PC8=scores_ht.scores[7],
        PC9=scores_ht.scores[8],
        PC10=scores_ht.scores[9])
    scatter(ht=scores_ht, x_location='PC1', y_location='PC2',
            plot_path=path_to_plots + '/SCATTER_1KG_PC1_PC2.png')
    scatter(ht=scores_ht, x_location='PC1', y_location='PC3',
            plot_path=path_to_plots + '/SCATTER_1KG_PC1_PC3.png')
    scatter(ht=scores_ht, x_location='PC2', y_location='PC3',
            plot_path=path_to_plots + '/SCATTER_1KG_PC2_PC3.png')
    scatter(ht=scores_ht, x_location='PC3', y_location='PC4',
            plot_path=path_to_plots + '/SCATTER_1KG_PC3_PC4.png')
    scatter(ht=scores_ht, x_location='PC4', y_location='PC5',
            plot_path=path_to_plots + '/SCATTER_1KG_PC4_PC5.png')
    scatter(ht=scores_ht, x_location='PC5', y_location='PC6',
            plot_path=path_to_plots + '/SCATTER_1KG_PC5_PC6.png')

    # Would be good to do some loadings plot to check whether the PCs captures batch effect or true structure


if __name__ == '__main__':
    path_to_mt = 'gs://unicorn-qc/post-imputation-qc/merged/mt/qc-af_0.005.mt'
    path_to_pruned_mt = 'gs://unicorn-qc/post-imputation-qc/pca/ld_pruned.mt'
    path_to_scores_ht = 'gs://unicorn-qc/post-imputation-qc/pca/scores.kt'
    path_to_loadings_ht = 'gs://unicorn-qc/post-imputation-qc/pca/loadings.kt'
    path_to_plots = 'gs://unicorn-qc/post-imputation-qc/summary'

    __main__(
        path_to_mt=path_to_mt,
        path_to_pruned_mt=path_to_pruned_mt,
        path_to_scores_ht=path_to_scores_ht,
        path_to_loadings_ht=path_to_loadings_ht,
        path_to_plots=path_to_plots)
