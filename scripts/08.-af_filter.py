import os
import sys
import hail as hl
from typing import *

arrs = [
	'info_Axiom_KP_UCSF_EUR',
	'info_Affymetrix_6.0',
	'info_Human550',
	'info_Human610660',
	'info_HumanOmni',
]

mt = hl.read_matrix_table('gs://unicorn-qc/post-imputation-qc/merged/mt/pre_qc.mt')
mt = mt.annotate_entries(MGP = hl.max(mt.GP))

# Filter to allele frequencies > 0.005
for arr in arrs:
	mt = mt.filter_rows(mt[arr].MAF > 0.005, keep = True)

print('There are {n_cols} samples '\
      'and {n_rows} variants '\
      'in post-imputation data'.format(
            n_cols = mt.count_cols(),
            n_rows = mt.count_rows()))
# There are 54125 samples and 8604906 variants in post-imputation data
mt.write('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-af_0.005.mt', overwrite=True)










