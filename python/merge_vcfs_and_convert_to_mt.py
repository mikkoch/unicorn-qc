import os
import sys
import hail as hl
from typing import *
from functools import reduce


def union_cols_all(mt_list) -> hl.MatrixTable:
    if type(mt_list) is list:
        mt = reduce(lambda x, y: x.union_cols(y), mt_list)
        return mt
    else:
        return mt_list


def merge_vcfs_and_convert_to_mt(path_to_vcfs: List[str], array_types: List[str], output_path: str):
     """
        This function is for merging vcf to MatrixTable. 
        ----
        :param List[str] path_to_vcfs: list of path to vcf files to be merged
        :param List[str] array_types: list of array types
        :param str output_path: output path for final matrixtable file
        """

    mt_list = list()
    variant_info = dict()
    variant_qc = dict()

    for path_to_vcf, array_type in zip(path_to_vcfs, array_types):

        mt = hl.import_vcf(path = path_to_vcf, min_partitions = 1000)
        mt = mt.annotate_cols(is_case = 0, array_type = array_type)
        mt = hl.variant_qc(mt)

        variant_info[array_type] =  mt.rows().select("info")
        variant_qc[array_type] =  mt.rows().select("variant_qc")

        mt_list.append(mt)

    mt = union_cols_all(mt_list) 

    for array_type in array_types:
        mt = mt.annotate_rows(**{'info_all.' + array_type: variant_info.get(array_type)[mt.locus, mt.alleles].info})
        mt = mt.annotate_rows(**{'qc_all.' + array_type: variant_qc.get(array_type)[mt.locus, mt.alleles].variant_qc})

    mt.write(output = output_path, overwrite = True)
