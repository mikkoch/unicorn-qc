import os
import sys
import hail as hl
from typing import *
from functools import reduce
sys.path.insert(0, '/home/danfengc/unicorn')


def __main__(path_to_vcfs: List[str], array_types:List[str], output_path:str):

    mt_list = list()
    variant_info_list = list()

    for path_to_vcf, array_type in zip(path_to_vcfs, array_types):

        mt = hl.import_vcf(path=path_to_vcf, min_partitions=1000)
        mt = mt.annotate_cols(
            is_case=0,
            array_type=array_type)

        # Row
        # fields:
        # 'locus': locus < GRCh37 >
        # 'alleles': array < str >
        # 'rsid': str
        # 'qual': float64
        # 'filters': set < str >
        # 'info': struct
        # {
        #     AF: array < float64 >,
        #     MAF: float64,
        #     R2: float64,
        #     ER2: float64,
        #     IMPUTED: bool,
        #     TYPED: bool,
        #     TYPED_ONLY: bool
        # }

        variant_info = mt.rows()
        variant_info = variant_info.transmute(**{
            "info_{array_type}".format(array_type=array_type): variant_info.info})

        mt_list.append(mt.select_rows())
        variant_info_list.append(variant_info)
    mt = reduce(lambda x, y: x.union_cols(y), mt_list)

    for variant_info in variant_info_list:
        mt = mt.annotate_rows(**variant_info[mt.row_key])

    mt.write(output=output_path, overwrite=True)

if __name__ == '__main__':
    array_types = [
        'Axiom_KP_UCSF_EUR',
        'Affymetrix_6.0',
        'Human550',
        'Human610660',
        'HumanOmni']
    path_to_vcfs = [
        'gs://unicorn-qc/post-imputation-qc/Axiom_KP_UCSF_EUR/vcf/chr*.dose.vcf.gz.bgz',
        'gs://unicorn-qc/post-imputation-qc/Affymetrix_6.0/vcf/chr*.dose.vcf.gz.bgz',
        'gs://unicorn-qc/post-imputation-qc/Human550/vcf/chr*.dose.vcf.gz.bgz',
        'gs://unicorn-qc/post-imputation-qc/Human610660/vcf/chr*.dose.vcf.gz.bgz',
        'gs://unicorn-qc/post-imputation-qc/HumanOmni/vcf/chr*.dose.vcf.gz.bgz']
    output_path = 'gs://unicorn-qc/post-imputation-qc/merged/mt/pre_qc.mt'

    __main__(
        path_to_vcfs=path_to_vcfs,
        array_types=array_types,
        output_path=output_path)