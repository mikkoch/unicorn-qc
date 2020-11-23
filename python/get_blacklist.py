import sys
import json

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *

def get_blacklist(path_to_mt: str, output_root_directory: path, array_types: List[str], 
	er2_thresh: float = 0.6):
	"""
        This function generates a blacklist of SNPs based on ER2
        ----
        :param str path_to_mt path to matrix table
        :param str output_root_directory
        :param List[str] list of string type 
        :param er2_thresh threshold for ER2
        """

    mt = hl.read_matrix_table(path_to_mt)
    directory_structure = build_cohort_qc_directory_structure(output_root_directory)

    for array_type in array_types:
    	mt = mt.filter_rows(mt['info_all.{0}'format(array_type)].ER2 <= er2_thresh)
        mt.rows().write(directory_structure["share"]+"/"+array_type+"-blacklist.tsv", overwrite = True)

   
