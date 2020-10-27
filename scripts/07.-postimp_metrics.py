import os
import sys
import numpy as np
from typing import *
sys.path.insert(0, '/home/danfengc/unicorn')
from qc.plotting import *

def __main__(
        path_to_mt: str,
        array_types:List[str],
        output_root:str,
        maximum_gp_thresh: float = 0.9,
        missingness_thresh: float = 0.1,
        missingness_thresh_for_pca: float = 0.01,
        R2_thresh_for_pca: float = 0.95,
        af_thresh_for_pca: float = 0.05,
        n_partitions_for_pruning: int=100):

    # Read MatrixTable from Google bucket
    mt = hl.read_matrix_table(path_to_mt)

    # show the basic information of the MatrixTable
    print(mt.describe())

    # posterior distribution 
    print('There are {n_cols} samples '\
          'and {n_rows} variants '\
          'in post-imputation data'.format(
            n_cols = mt.count_cols(),
            n_rows = mt.count_rows()))
    # There are 54125 samples and 39117105 variants in post-imputation data

    # show the entry fields
    print(mt.entries().describe())
    # Entry fields:
    # 'GT': call
    # 'DS': float64
    # 'HDS': array<float64>
    # 'GP': array<float64>

    # annotate maximum posterior probability
    mt = mt.annotate_entries(
        MGP = hl.max(mt.GP))

    # sample 
    small_mt = mt.sample_rows(0.01)

    # calculate the summary statistics for GP
    mgp_summ = small_mt.aggregate_entries(hl.agg.stats(small_mt.MGP))
    print(mgp_summ)

    # restrict to AF > 0.01, the distribution of maximum posterior probability
    for arr in array_types:
        small_mt = small_mt.filter_rows(
            small_mt['info_'+arr].MAF > 0.01)

    # calculate the summary statistics for GP
    mgp_summ = small_mt.aggregate_entries(hl.agg.stats(small_mt.MGP))
    print(mgp_summ)

    for thresh in np.linspace(0.3, 1, 8):
        print(small_mt.aggregate_entries(hl.agg.fraction(small_mt.MGP >= thresh)))

    # annotate missing 
    mt = mt.annotate_entries(
        is_missing = mt.MPG < maximum_gp_thresh, 
        is_called  = mt.MPG >=maximum_gp_thresh)
    mt = mt.annotate_rows(
        missing_rate = hl.agg.fraction(mt.is_missing),
        call_rate = hl.agg.fraction(mt.is_called))

    # calculate the distribution 



    # draw the distribution of minor allele frequencies & INFO scores 

    # draw the distribution of max posterior probability 

    # I should 

if __name__ == '__main__':
    array_types = [
        'Affymetrix_6.0',
        'Human550',
        'Human610660',
        'HumanOmni',
        'Axiom_KP_UCSF_EUR']
    path_to_mt  = 'gs://unicorn-qc/post-imputation-qc/merged/mt/pre_qc.mt'
    output_root = 'gs://unicorn-qc/post-imputation-qc/merged'

    __main__(
        path_to_mt=path_to_mt,
        array_types=array_types,
        output_root=output_root)
