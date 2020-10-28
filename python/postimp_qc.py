import sys
import json

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *

def postimp_qc(path_to_mt: str, output_root_directory: path, array_types: List[str], 
	info_score_thresh: float, hwe_thresh: float, maf_thresh: float):
	"""
        This function performs post-imputation QC including the following QC filters:
        - minor allele frequencies 
        - hardy weinberg equilibrium
        - imputation info score 
        ----
        :param str path_to_mt path to matrix table
        :param str output_root_directory
        :param List[str] list of string type 
        :param float info_score_thresh threshold for imputation info score
        :param hwe_thresh hardy weinberg equilibrium threshold 
        """

    mt = hl.read_matrix_table(path_to_mt)
    directory_structure = build_cohort_qc_directory_structure(output_root_directory)
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # Perform sample QC
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)

    for array_type in array_types:
    	mt = mt.filter_rows(mt['info_all.{0}'format(array_type)].MAF > maf_thresh)
    	mt = mt.filter_rows(mt['info_all.{0}'format(array_type)].R2  > info_score_thresh)
    	mt = mt.filter_rows(mt['qc_all.{0}'format(array_type)].het_freq_hwe > hwe_thresh)

    mt = mt.filter_rows((mt.variant_qc.AF[0] > maf_thresh) &
                        (mt.variant_qc.AF[1] > maf_thresh) &
                        (mt.variant_qc.het_freq_hwe > hwe_thresh))

    n_variants_after_qc, n_samples_after_qc = mt.count()
    # create a dictionary storing the meta information
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc}

    mt.write(directory_structure["share"]+"/all.mt", overwrite = True)

    # Export meta
    with hl.hadoop_open(directory_structure["summary"] + "/post_imputation_qc_summary.json", 'w') as outfile:
        json.dump(meta, outfile)
