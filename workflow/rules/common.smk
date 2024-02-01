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
def metadata_builder(file, SRS_dict = {}, discrepancy = False):
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
					    	  'paired':True if info[2]=='PAIRED' else False, 
					          'organism':info[3].replace(' ', '_'),
		            	      'tech':info[4],
						      'UMI':True if info[5].upper()=='YES' else False,
							  'Study': info[6]}
			else:
				# this is mostly for SRA having the 'paired' status wrong
				# don't want to hand-edit the main metadata file
				# so I think better to create a new file with
				# hand edited values for just the ones i want to change
				if discrepancy:
					runs = SRS_dict[SRS]['SRR']
					SRS_dict[SRS] = {'SRR':runs,
									 'paired':True if info[2]=='PAIRED' else False,
									 'organism':info[3],
									 'tech':info[4],
									 'UMI':True if info[5]=='YES' else False,
									 'Study': info[6]}
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


def get_whitelist_from_tech(tech):
	if tech == 'DropSeq':
		return f'references/whitelist/DropSeq/barcodes.txt'
	else:
		return('references/whitelist/10x/{tech}.txt')

def get_kallisto_quant_layout_flag(is_paired):
	if is_paired:
		return ''
	else: 
		return '--single -l 200 -s 30'

def DROPSEQ_samples_from_reference(quant_path ,srs_dict):
	out= []
	for sample in srs_dict.keys():
		if srs_dict[sample]['tech'] == 'DropSeq' :
			if srs_dict[sample]['organism'] == 'Homo_sapiens':
				out.append(f'{quant_path}/quant/{sample}/DropSeq/hs-homo_sapiens/output.sorted.bus')
			else:
				out.append(f'{quant_path}/quant/{sample}/DropSeq/mm-mus_musculus/output.sorted.bus')
	return out 

