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


def convert_to_vcf(mt:hl.MatrixTable, basename:str, output_root_for_unzip_vcf:str, output_root_for_zip_vcf:str,
                   temp_root: str='file:///home/danfengc'):
    """
        This function is for converting MatrixTable to vcf. It requires bgzip & vcf-sort to be installed
        through the following commands
            ```
            sudo apt-get install vcftools -y
            sudo apt-get install tabix -y
            ```
        ----
        :param MatrixTable mt: input hail matrix table
        :param str basename: basename of output vcf
        :param str output_root_for_unzip_vcf: output vcf
        :param str tmp_root: directory in master node for saving the temporary vcf files
        """
    for i in range(1, 23):
        dataset = mt.filter_rows(mt.locus.contig == hl.literal(str(i)))
        out = create_vcfname(output_root=output_root_for_unzip_vcf, basename="final", chromosome=i)
        hl.export_vcf(dataset, out)

        src = create_vcfname(output_root=output_root_for_unzip_vcf, basename="final", chromosome=i)
        des = create_vcfname(output_root=temp_root, basename="final", chromosome=i)
        hl.utils.hadoop_copy(src, des)

        command = "sudo /usr/bin/vcf-sort {0} | /usr/bin/bgzip -c > {1}".format(
            create_vcfname(output_root=temp_root.replace('file://', ''), basename=basename, chromosome=i),
            create_zipped_vcfname(output_root=temp_root.replace('file://', ''), basename=basename, chromosome=i))
        subprocess.call(command, shell=True)

        src = create_zipped_vcfname(output_root=temp_root, basename=basename, chromosome=i)
        des = create_zipped_vcfname(output_root=output_root_for_zip_vcf, basename=basename, chromosome=i)
        hl.utils.hadoop_copy(src, des)