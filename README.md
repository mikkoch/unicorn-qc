## A Data Harmonization Pipeline To Leverage External Controls AND BOOST POWER IN GWAS
 

We propose a unified pipeline to harmonize control samples from different cohort studies that may have been genotyped using multiple different array platforms. The pipeline contains four modules: (i) *Stratification* Within-array processing, (ii) *Imputation*, (iii) *Cross-array comparison*, and (iv) *Re-imputation*.

### Module 1: Stratification. 

- *Cohort-level QC:* perform cohort-level quality control [cohort_qc.py](python/cohort_qc.py)
- *Ancestry matching:* infer the ancestry of each sample [match_ancestry.py](python/match_ancestry.py)
- *Merging:* merge samples sharing the same genotyping array and ancestry group [merge_cohorts.py](python/merge_cohorts.py)
- *Array-level QC:* perform array-level quality control  [array_qc.py](python/array_qc.py)

### Module 2: Imputation. 

We use the [Michigan Imputation server](https://imputationserver.sph.umich.edu/index.html), with 1000 Genomes data as the reference panel. It requires vcf-formatted input, which could be obtained with the script [convert\_to_vcf.py](python/convert_to_vcf.py)
.
 
### Module 3: Cross-array comparison. 
- *Merging:* convert and merge imputed VCFs from Michigan imputation server to MatrixTables [merge\_vcfs_and\_convert\_to\_mt.py](python/merge_vcfs_and_convert_to_mt.py)
- *Post-imputation QC.:* removes variants with low minor allele frequencies, or small Hardy-Weinberg Equilibrium p-values, or low imputation info scores [postimp_qc.py](python/postimp_qc.py)
- *Cross-array pseudo GWAS*: perform cross-array type pseudo-GWAS [cross\_array_comparison.py](python/cross_array_comparison.py)


### Module 4: Re-imputation. 
- *Generating Blacklist of SNPs:* generate a blacklist of SNPs based on ER2 filter [get_blacklist.py](python/get_blacklist.py)
- *Re-run analysis.* remove the blacklist of SNPs from the data set produced by module 1, and then re-run module 2, 3, 4.
 
<br/><br/>

### Blacklist 
- We applied this pipeline to aggregate 27,517 European samples from 16 collections within dbGaP. Here we provide the list of problematic SNPs identified and removed in different steps of the pipeline.  

- The blacklist of SNPs based on ER2 filter for array types [HumanHap300](blacklist/Human300.tsv), [HumanHap550](blacklist/Human550.tsv), [HumanHap610](blacklist/Human610.tsv), [HumanHap660](blacklist/Human660.tsv), [Affymetrix 6.0](blacklist/Affymetrix_6.tsv) and [Axiom_KP_UCSF_EUR](blacklist/Axiom_KP_UCSF_EUR.tsv)
- The blacklist of SNPs that exhibit significant p-values in cross-array-type comparison: [cross-comparison blacklist](blacklist/cross_comparison.tsv)

<br/><br/>
 
Reference: 

Chen, Tashman, Palmer, Neale, Roeder, Bloemendal, Churchhouse and Ke (2020). A data harmonization pipeline to leverage external controls and boost power in GWAS. 




