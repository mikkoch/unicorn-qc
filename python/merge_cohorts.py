import sys
import json
from typing import *
from functools import reduce

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *


def merge_cohorts(mt_path_list: List[str], cohort_labels: List[str],
    population: str, array: str, output_root_directory: str):
    """
        This function merges samples sharing the same genotyping array and ancestry group.
        ----
        :param List[str] mt_path_list: a list of path to MatrixTable from same genotyping arrays
        :param List[str] cohort_labels: a list of cohort labels
        :param str output_root_directory
        :param str population
        :param str array
        """

    # Merge MatrixTable from same population and same genotyping array
    mt_list = list()
    for mt_path, cohort in zip(mt_path_list, cohort_labels):
        mt = hl.read_matrix_table(mt_path)
        mt = mt.annotate_cols(cohort=cohort)
        mt_list.append(mt)
    mt = reduce(lambda x, y: x.union_cols(y), mt_list)
    mt = mt.annotate_cols(population = population, array = array)

    # Build Directory structure
    directory_structure = build_array_qc_directory_structure(output_root_directory, population=population, array=array)

    # Write QC'ed data to google bucket
    mt.write(directory_structure["share"]+"/final.mt", overwrite = True)
    