## A Data Harmonization Pipeline To Leverage External Controls 

We propose a unified pipeline to harmonize control samples from different cohort studies that may have been genotyped using multiple different array platforms. The pipeline contains four modules: (i) *Stratification*, (ii) *Imputation*, (iii) *Cross-array comparison*, and (iv) *Re-imputation*.

### Module 1: Stratification. 

- *Cohort-level QC:* perform cohort-level quality control: [cohort_qc.py](python/cohort_qc.py)
- *Ancestry matching:* infer the ancestry of each sample: [match_ancestry.py](python/match_ancestry.py)
- *Merging:* merge samples sharing the same genotyping array and ancestry group [merge_cohorts.py](python/merge_cohorts.py)
- *Array-level QC:* perform array-level quality control:  [array_qc.py](python/array_qc.py)

### Module 2: Imputation. 

We use the [Michigan Imputation server](https://imputationserver.sph.umich.edu/index.html), with 1000 Genomes data as the reference panel. It requires vcf-formatted input, which could be obtained with the script [convert\_to_vcf.py](python/convert_to_vcf.py)
.
 
### Module 3: Cross-array comparison.

### Module 4: Re-imputation. 
 



