## UNICORN QC PIPELINE

### STEP 0: ALIGNING DATA TO HRC 

### STEP 1: POPULATION REFERENCE 

This step is to set up the population reference (1000 genomes and Ashkenazi Jewish samples) for ancestry assignment. 

**1000 Genomes** 

We use the QC'ed 1000 genomes for population assignment (provided by Stephan). The reasons for doing that are: 

1. projecting onto raw 1000 genomes leads to severe shrinkage issue 
2. PC1 calculated from joint PCA with raw 1000 genomes captures batch effects

```
## Broad Server
/psych/genetics_data/ripke/references_outdated/hapmap_ref/impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/4pops/qc/pop_4pop_mix_SEQ.*
/psych/genetics_data/ripke/references_outdated/hapmap_ref/impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/4pops/qc/pop_euro_eur_SEQ.*
## Google Bucket
gs://unicorn-resources/1000_genomes/pop_4pop_mix_SEQ.*
gs://unicorn-resources/1000_genomes/pop_euro_eur_SEQ.*
```

**Ashkenazi Jewish Genotypes**

We use controls sample from a Ashkenazi GWAS cohort as reference population.

```
gs://unicorn-resources/Ashkenazi_Jewish_Samples/Ashkenazi_Jewish_Samples.*

```

**Steps**

To convert the population reference in PLINK format to hail MatrixTable that we can apply to downstream QC, the following steps are required: 

1. Align population references to HRC
2. Convert population references in PLINK format to Hail MatrixTable
3. Merge 1KG EUR samples with AJ samples 

### STEP 2: COHORT QC 

This step is to perform QC at cohort level; the QC metrics are listed below. This step also performs ancestry inference, and assign each sample to one of the `AMR`, `AFR`, `EAS`, `EUR(mainland)`, `AJ` and `FIN` populations. The input is genotype in PLINK format that has been aligned to HRC reference. Tht output is path to the directory for saving the QC results.

**QC metrics**

```
1. Sample call rate > 0.98
2. Per-chromosome sample call rate > 0.50
3. Inbreeding coefficient should be within 3 standard deviation from its mean 
4. Remove related samples (pi-hat threshold = 0.0625)
5. Remove sex-error samples
6. Variant call rate > 0.99
7. Minor allele frequencies > 0.01
8. Within each population, p-hwe > 1e-4
9. Autosome only
10. Skip any genotyping file that sample size is less than 50
```

**Ancestry Inference**

We calculate PCA using population references (1000 genomes), and then project genotypes onto the PC space. We then uses a random forest model to assign population labels based on the results of PCA. 

**Output**

```
1. qc'ed hail matrix table 
2. histograms of qc metrics
3. PCA plot for ancestry assignment 
4. json file that summarizes QC process
```
Please check the following google bucket as an example output `gs://unicorn-qc/pre-imputation-qc/cohort-qc/geneva_t2d_hpfs`

### STEP 3: ARRAY QC

This step is to perform QC at array level including: 

1. Merge cohorts genotyped with the same array and from same population
2. Perform QC at the array level
3. convert the QC'ed hail matrix table to VCF

**QC metrics**

```
1. Sample call rate > 0.98
2. Per-chromosome sample call rate > 0.50
3. Remove related samples (pi-hat threshold = 0.0625)
4. Variant call rate > 0.99
5. Minor allele frequencies > 0.01
6. Within each population, p-hwe > 1e-4
7. Pseudo case-control GWAS (tagging samples from one cohort as cases, from the other cohorts as controls) p-values > 1e-4 
8. Pseudo case-control GWAS against 1KG p-values > 1e-4
```

**Output**

```
1. VCF 
2. histograms of qc metrics
3. PCA plot for ancestry assignment 
4. json file that summarizes QC process
```

Please check the following google bucket as an example output `unicorn-qc/pre-imputation-qc/array-qc/Affymetrix6`

### HOW TO RUN

First of all, we set up the dataproc using `cloudtools` specifying the additional python modules `scipy` and `scikit-learn`. For instance

```
cluster start tst \
    --version 0.2 \
    --spark 2.2.0 \
    --packages scikit-learn,sympy
```

Copy scripts to the master node by

```
gcloud compute ssh tst-m 
mkdir unicorn
logout
gcloud compute scp --recurse ${PATH_TO_UNICORN}/main tst-m:/home/danfengc/unicorn/main
gcloud compute scp --recurse ${PATH_TO_UNICORN}/qc tst-m:/home/danfengc/unicorn/qc
```

**Step 0: Align data to HRC**

The google bucket contains all raw control samples: 

```
gs://control_repo_by_platform
```

The PLINK bfiles that need to be proceed are the following: 

```
gs://control_repo_by_platform/Affymetrix_6.0/geneva_t2d/jan_24_2018_geneva_t2d_hpfs_zeroed_out.updated.*
gs://control_repo_by_platform/Affymetrix_6.0/geneva_t2d/jan_24_2018_geneva_t2d_nhs_zeroed_out.updated.*
gs://control_repo_by_platform/Affymetrix_6.0/gwas_of_schizophrenia/dec_14_2017_gwas_of_schizophrenia_european_filtered.*
gs://control_repo_by_platform/Affymetrix_6.0/neurodevelopmental_genomics/jan_22_2018_neuro_develop_affy6.0.*
gs://control_repo_by_platform/Axiom_Genome-Wide_Human_Origins/neurodevelopmental_genomics/jan_22_2018_neuro_develop_axiom.*
gs://control_repo_by_platform/Axiom_Genome-Wide_Human_Origins_Tx/neurodevelopmental_genomics/jan_22_2018_neuro_develop_axiom_Tx.*
gs://control_repo_by_platform/Axiom_KP_UCSF_EUR/gera/eur_GERA.*
gs://control_repo_by_platform/Human1M-Duov3_B/neurodevelopmental_genomics/jan_22_2018_neuro_develop_human1m.*
gs://control_repo_by_platform/Human610W-Quad_v1_B/ad_family_study/dec_18_2017_ad_family_study.*
gs://control_repo_by_platform/Human610W-Quad_v1_B/neurodevelopmental_genomics/jan_22_2018_neuro_develop_human610_quadv1_version1.*
gs://control_repo_by_platform/Human610W-Quad_v1_B/neurodevelopmental_genomics/jan_22_2018_neuro_develop_human610_quadv1_version2.*
gs://control_repo_by_platform/Human610W-Quad_v1_B/panscan/dec_18_2017_panscan_human610_quadv1_c1.*
gs://control_repo_by_platform/Human610W-Quad_v1_B/study_of_pediatric_disorders/jan_19_2018_controls_610_quadv1.*
gs://control_repo_by_platform/Human660W-Quad_v1_A/glaucoma_gwas/jan_18_2018_glaucoma_gwas_subject_level_filtered.*
gs://control_repo_by_platform/Human660W-Quad_v1_A/vte_gwas/jan_18_2018_vts_gwas_subject_level.*
gs://control_repo_by_platform/HumanCNV370_v1/age_related_macular_degeneration/jan_22_2018_amd_mmap_cohort.*
gs://control_repo_by_platform/HumanHap550_v1.1/exome_sequencing_of_als/jan_17_2018_als_exome_seq_v1_controls_550.*
gs://control_repo_by_platform/HumanHap550_v1.1/exome_sequencing_of_als/jan_17_2018_als_exome_seq_v2_italian_550v1.*
gs://control_repo_by_platform/HumanHap550_v1.1/neurodevelopmental_genomics/jan_22_2018_neuro_develop_humanhap550_v1.*
gs://control_repo_by_platform/HumanHap550_v1.1/study_of_pediatric_disorders/jan_19_2018_controls_bdchp_1x10_550.*
gs://control_repo_by_platform/HumanHap550_v3.0/exome_sequencing_of_als/jan_17_2018_als_exome_seq_v2_controls_us_550v3.*
gs://control_repo_by_platform/HumanHap550_v3.0/neurodevelopmental_genomics/jan_22_2018_neuro_develop_human550_v3_version1.*
gs://control_repo_by_platform/HumanHap550_v3.0/neurodevelopmental_genomics/jan_22_2018_neuro_develop_human550_v3_version2.*
gs://control_repo_by_platform/HumanHap550_v3.0/panscan/jan_24_2018_panscan_human550v3_c1.case.control.*
gs://control_repo_by_platform/HumanHap550_v3.0/panscan/jan_24_2018_panscan_human550v3_c1.cohort.*
gs://control_repo_by_platform/HumanHap550_v3.0/study_of_pediatric_disorders/jan_19_2018_controls_550v3.*
gs://control_repo_by_platform/HumanOmni5-Quad/stroke_genetics_network/jan_22_2018_stroke_network_subject_level.*
gs://control_repo_by_platform/HumanOmniExpress_12v1.0A/neurodevelopmental_genomics/jan_22_2018_humanOmniExpress_12v1a_version1.*
gs://control_repo_by_platform/HumanOmniExpress_12v1.0A/neurodevelopmental_genomics/jan_22_2018_humanOmniExpress_12v1a_version2.*
gs://control_repo_by_platform/HumanOmniExpress_12v1.0A/panscan/dec_18_2017_panscan_human_omni_express12.*
gs://control_repo_by_platform/HumanOmniExpress_12v1.1B/neurodevelopmental_genomics/jan_22_2018_humanOmniExpress_12v1b.* 
```

**Step 1: population reference**

The population reference files that need to be aligned to HRC are here: 

```
gs://unicorn-resources/1000_genomes/pop_4pop_mix_SEQ.*
gs://unicorn-resources/1000_genomes/pop_euro_eur_SEQ.*
gs://unicorn-resources/Ashkenazi_Jewish_Samples/Ashkenazi_Jewish_Samples.*
```

And then convert them to MatrixTable for downstream analysis. 

