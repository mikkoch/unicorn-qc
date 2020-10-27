import sys
import json

from qc.utils import *
from qc.management import *
from qc.plotting import *
from main.pca import *

def cohort_qc(bfile: str, output_root_directory: str, sample_call_rate_thresh: float, per_chr_sample_call_rate_thresh: float,
	variant_call_rate_thresh: float, maf_thresh: float,  n_partitions: int = 200):
    """
        This function performs cohort-level QC including the following QC filters:
        - sample call rate
        - sample call rate for each chromosome 
        - variant call rate 
        - inbreeding coefficient should be less than 3 std from its mean
        - remove related samples (pi-hcat threshold = 0.0625)
        - remove sex-error samples
        - minor allele frequencies 
        - autosome only
        ----
        :param str bfile: path to input PLINK bfile
        :param str output_root_directory
        :param float sample_call_rate_thresh: sample call rate threshold
        :param float per_chr_sample_call_rate_thresh: per chromosome sample call rate
        :param float variant_call_rate_thresh: variant call rate threshold
        :param float maf_thresh: minor allele frequencies threshold
        :param n_partitions: number of partitions
        """

    mt = hl.import_plink(bed=bfile+'.bed', bim=bfile+'.bim', fam=bfile+'.fam', min_partitions=n_partitions)
    directory_structure = build_cohort_qc_directory_structure(output_root_directory)
    n_variants_before_qc, n_samples_before_qc = mt.count()

    # Perform sample QC
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)
    mt = calculate_per_chr_sample_call_rate(mt = mt, sample_call_rate_col = 'scr_chr')  # calculate per-chromosome call rate for each sample
    mt = calculate_inbreeding_coefficients(mt = mt, ib_col = 'ib', ib_global = 'ib_stats')  # calculate inbreeding coefficient
    mt = mt.cache()

    # Calculate number of samples & variants that fail each filter
    n_sample_call_rate_thresh = mt.aggregate_cols(hl.agg.count_where(mt.sample_qc.call_rate < sample_call_rate_thresh))
    n_per_chr_sample_call_rate_thresh = mt.aggregate_cols(hl.agg.count_where(mt.scr_chr < per_chr_sample_call_rate_thresh))
    n_inbreeding = mt.aggregate_cols(hl.agg.count_where(hl.abs(mt.ib.f_stat-mt.ib_stats.mean) > 3 * mt.ib_stats.stdev))
    n_variant_call_rate_thresh  = mt.aggregate_rows(hl.agg.count_where(mt.variant_qc.call_rate < variant_call_rate_thresh))
    n_maf_thresh = mt.aggregate_rows(hl.agg.count_where((mt.variant_qc.AF[0] < maf_thresh) | (mt.variant_qc.AF[1] < maf_thresh)))

    # Filter samples & variants
    mt = mt.filter_cols((mt.sample_qc.call_rate > sample_call_rate_thresh) &
                        (hl.abs(mt.ib.f_stat - mt.ib_stats.mean) <= 3 * mt.ib_stats.stdev) &
                        (mt.scr_chr > per_chr_sample_call_rate_thresh), keep=True)
    mt = mt.filter_rows((mt.variant_qc.call_rate > variant_call_rate_thresh) &
                        (mt.variant_qc.AF[0] > maf_thresh) &
                        (mt.variant_qc.AF[1] > maf_thresh) &
                        (mt.locus.in_autosome()), keep=True)

    # Count samples & variants after QC
    n_variants_after_qc, n_samples_after_qc = mt.count()

    # aggregate summary statistics
    meta = {"Number of Samples before QC": n_samples_before_qc,
            "Number of Samples after QC:": n_samples_after_qc,
            "Number of Variants before QC": n_variants_before_qc,
            "Number of Variants after QC": n_variants_after_qc,
            "Sample QC Summary": {
                "IDs: call rate < %s" % n_sample_call_rate_thresh: n_sample_call_rate_thresh,
                "IDs: minimum per-chromosome call rate < %s" %n_per_chr_sample_call_rate_thresh:n_per_chr_sample_call_rate_thresh,
                "IDs: inbreeding coefficient is 3 stdev from its mean": n_inbreeding},
            "Variant QC Summary": {
                "SNPs: call rate < %s" % n_variant_call_rate_thresh: n_variant_call_rate_thresh,
                "SNPs: minor allele frequency < %s" % n_maf_thresh: n_maf_thresh}
            }
    mt.write(directory_structure["share"]+"/all.mt", overwrite=True)

    # export meta
    with hl.hadoop_open(directory_structure["summary"] + "/cohort_qc_summary.json", 'w') as outfile:
        json.dump(meta, outfile)



