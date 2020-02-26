import os
import sys
import hail as hl
from typing import *


mt = hl.read_matrix_table('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-af_0.005.mt')


# print the summary statistics of maximum posterior probability
summary = mt.aggregate_entries(hl.agg.stats(mt.MGP))
print(summary)
# mean=0.9859765560738231, stdev=0.055279653337162775, min=0.444, max=1.0, n=4646631250, sum=4581469477.220003

for thresh in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
	print('There are {n} genotypes with '\
		  'maximum posterior probability '\
		  'lower than {thresh}'.format(
		  		n = mt.aggregate_entries(hl.agg.count_where(mt.MGP < thresh)), 
		  		thresh = thresh))


# convert to missing if maximum posterior probability lower than 0.9
mt = hl.read_matrix_table('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-af_0.005.mt')
mt = mt.transmute_entries(
	GT = hl.cond(mt.MGP < 0.9, hl.null(hl.tcall), mt.MGP))
mt.write('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-mgp_0.9.mt', overwrite=True)


# convert to missing if maximum posterior probability lower than 0.95
mt = hl.read_matrix_table('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-af_0.005.mt')
mt = mt.transmute_entries(
	GT = hl.cond(mt.MGP < 0.95, hl.null(hl.tcall), mt.MGP))
mt.write('gs://unicorn-qc/post-imputation-qc/merged/mt/qc-mgp_0.95.mt', overwrite=True)





