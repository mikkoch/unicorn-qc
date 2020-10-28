import os
import subprocess
import hail as hl


def build_cohort_qc_directory_structure(dir_name:str)->dict:
    """
        This function is to create the directory structure for saving the
        cohort QC results.
        ----
        :param str dir_name: output root directory
        """
    return {"main": dir_name,
            "share":os.path.join(dir_name, 'share'),
            "pca": os.path.join(dir_name, 'pca'),
            "summary":os.path.join(dir_name, 'summary')}


def build_array_qc_directory_structure(dir_name:str, population:str, array:str)->dict:
    """
        This function is to create the directory structure for saving the
        cohort QC results.
        ----
        :param str dir_name: output root directory
        :param str population
        :param str array
        """
    dir_name = os.path.join(dir_name, array, population)
    return {
        "main": dir_name,
        "share":os.path.join(dir_name, 'share'),
        "pca": os.path.join(dir_name, 'pca'),
        "summary": os.path.join(dir_name, 'summary'),
        "vcf": os.path.join(dir_name, 'vcf')}


def create_zipped_vcfname(output_root: str, basename: str, chromosome: int) -> str:
    """
        helper function for `convert_to_vcf`
        """
    zipped_vcfname = os.path.join(output_root, basename + ".chr{0}.vcf.gz".format(chromosome))
    return zipped_vcfname


def create_vcfname(output_root: str, basename: str, chromosome: int) -> str:
    """
        helper function for `convert_to_vcf`
        """
    vcfname = os.path.join(output_root, basename + ".chr{0}.vcf".format(chromosome))
    return vcfname

