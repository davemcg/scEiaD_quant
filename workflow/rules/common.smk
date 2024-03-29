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
import gzip

def opener(filename, mode):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode)
	else:
		return open(filename, mode)
# builds dictionary of dictionaries where first dict key is SRS 
# and second dict key are important SRS properties
def sample_dict_builder(filename, SRS_dict = {}):
	with opener(filename, 'rt') as txt:
		for line in txt:
			if line[0] == '#':
				continue
			info = line.strip('\n').split('\t')
			if info[0] == 'sample_accession':
				continue
			SRS = info[0]
			if SRS not in SRS_dict:
				SRS_dict[SRS]={'SRR': [info[1]],
							  'paired': True if info[2] == 'PAIRED' else False,
							  'parity':' --parity paired ' if info[2]=='PAIRED' else '--parity single --fragment-l 200 --fragment-s 30 --tcc ',
							  'tech':info[4]}
			else:
				runs = SRS_dict[SRS]['SRR']
				runs.append(info[1])
				SRS_dict[SRS]['SRR'] = runs
	return(SRS_dict)

def make_core_outputs(filename, quant_path, umi_outputs = [], bulk_outputs = []):
	with opener(filename, 'rt') as txt:
		for line in txt:
			if line[0] == '#':
				continue
			info = line.strip('\n').split('\t')
			if info[0] == 'sample_accession':
				continue
			if info[5].upper() == 'TRUE':
				umi_outputs.append(quant_path + '/quant/' + info[0] + \
							'/' + info[3] + '/' + info[6] + '_' + info[7])
			else:
				bulk_outputs.append(quant_path + '/quant/' + info[0] + \
							'/' + info[3] + '/' + info[6] + '_' + info[7])
				
		return([list(set(umi_outputs)), list(set(bulk_outputs))])

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

git_dir = config['git_dir']

#wildcard_constraints:
#	workflow = ['nac|standard'],
#	sum = ['total|nucleus|cell']
