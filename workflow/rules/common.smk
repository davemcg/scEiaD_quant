print("Loading common.smk")
import pprint
pp = pprint.PrettyPrinter(width=41, compact=True) 
import subprocess as sp
import tempfile
import yaml
import json
import string
import itertools
import time 
import re

# builds dictionary of dictionaries where first dict key is SRS 
# and second dict key are SRS properties
def metadata_builder(file, SRS_dict = {}):
	with open(file) as file:
		for line in file:
			if line[0] == '#':
				continue
			info = line.strip('\n').split('\t')
			if info[0] == 'sample_accession':
				continue
			SRS = info[0]
			if SRS not in SRS_dict:
				SRS_dict[SRS]={'SRR': [info[1]],
							  'paired': True if info[2] == 'PAIRED' else False,
							  'parity':' --parity paired ' if info[2]=='PAIRED' else '--parity single --fragment-l 200 --fragment-s 30 ',
							  'ref':info[3],
							  'tech':info[4],
							  'umi': True if info[5] == 'TRUE' else False,
							  'workflow': info[6]}
			else:
				runs = SRS_dict[SRS]['SRR']
				runs.append(info[1])
				SRS_dict[SRS]['SRR'] = runs
	return(SRS_dict)

def lookup_run_from_SRS(SRS, fqp):
	SRR_files=SRS_dict[SRS]['SRR']
	out = []
	for SRR in SRR_files:
		if SRS_dict[SRS]['paired']:
			#PE
			out.append(f'{fqp}/fastq/{SRR}_1.fastq.gz')
			out.append(f'{fqp}/fastq/{SRR}_2.fastq.gz')
		else:
			#SE
			out.append(f'{fqp}/fastq/{SRR}.fastq.gz')
	return(out)

def ref_builder(file, ref_dict = {}):
	with open(file) as file:
		for line in file:
			info = line.strip('\n').split('\t')
			ref = info[0]
			ref_dict[ref]={'gtf': [info[1]],
							'fasta': info[2]}
	return(ref_dict)
