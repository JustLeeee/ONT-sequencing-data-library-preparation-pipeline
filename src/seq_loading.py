import argparse
import numpy as np
import re
import random
import time
import sys
import datetime
from Bio import SeqIO

############################################################################################################################
#------------- ----------------------------functions ---------------------------------------------------#
'''
Loading raw sequence files in fasta or txt format. Each sequence will generate a fasta format file.
'''
def load_genome(input_file):
	'''
	Note: There must be no newline characters in each sequence.
	'''
	id_list = list()
	sequence_list = list()
	for record in SeqIO.parse(input_file, "fasta"):
		id_list.append(record.id)
		sequence_list.append(record.seq)

	sequence = list(map(str, sequence_list))
	with open(input_file, 'r') as f:
		text = f.read()
		lines = text.splitlines()
	headers = filter(lambda x: '>' in x, lines)
	headers = list(headers)
	seq_lens = list(map(len, sequence))
	sequence = ''.join(sequence)
	return sequence, headers, seq_lens

def save_genome(genome, output_file, headers, seq_lens):
	'''
	Result: Each sequence will generate a  file in fasta format.
	The content of each file is as follows:
	'Seq Index': >***************..........
	'Sequence' : ATCGATCGATCG..............
	'''
	seq_lens = [0]+seq_lens
	for i in range(len(headers)):
		output_file1 = output_file + '/load_processing_seq_{}'.format(i) + '.fasta'
		with open(output_file1, 'w') as f:
			f.write(headers[i]+'\n')
			sequence = genome[sum(seq_lens[:i+1]):sum(seq_lens[:i+2])]
			f.write(sequence+'\n')

############################################################################################################################
# ----------------------------------------------- main -----------------------------------------------#
def initialization_parameters():
	parser = argparse.ArgumentParser()
	parser.add_argument('-input', type=str, required=True, help='The input file path')
	#parser.add_argument('-seq_len', type=int, required=True,help='Specifies the length of the sequence.')
	parser.add_argument('-output', type=str, required=True,help='The output file for storing qualified sequences.')
	args = parser.parse_args()
	return args


def main():
	args = initialization_parameters()
	a,b,c = load_genome(args.input)
	save_genome(a,args.output,b,c)
	#print('End')

##############################################################################################################################
if __name__ == "__main__":
	time_start = time.time()

	main()
	# print('End')

	time_end = time.time()
	time_sum = time_end - time_start

	print('Program running time: %d' % time_sum)