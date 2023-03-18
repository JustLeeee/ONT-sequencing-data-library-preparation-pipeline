import argparse
import numpy as np
import scipy.stats as st
import random
import time
import sys
import datetime

#########################################################################################################################
#---------------------------------The length distribution for the reads sampling--------------------------------#
#----------Part1: Three read length distribution models-------------#

'''
A total of three distribution models are set up in this part. 
Users can choose one model according to their own needs of  their raw datas. 
Users can also use their own defined distribution model through the interface in these codes.
'''

# s18
# beta distribution with parameters
# 1.7780912583853339,
# 7.8927794827841975,
# 316.75817723758286,
# 34191.25716704056
def draw_beta_dis(size, mean, seed):
	'''
	mean: the rough mean of the read length, default=8000
	return: 'samples' is a list.
			The elements in the list are the length of each reading to be generated
	'''
	samples = st.beta.rvs(1.778, 7.892, 316.758,
		34191.257, size=size, random_state=seed)
	samples = samples*mean/6615.0
	samples = samples.astype(int)
	samples = np.clip(samples, 1, len(genome))
	return samples

# human
# alpha distribution with parameters
# 0.0058193182047962533,
# -49.180482198937398,
# 1663.9103931473874
def draw_alpha_dis(size, mean, seed):
	'''
	mean: the rough mean of the read length, default=8000
	return: 'samples' is a list.
			The elements in the list are the length of each reading to be generated
	'''
	samples = st.alpha.rvs(0.00582,-49.1805,1663.91,
		size=size, random_state=seed)
	samples = samples*mean/7106.0
	samples = samples.astype(int)
	samples = np.clip(samples, 1, len(genome))
	return samples


def draw_expon_dis(size, mean, seed):
	'''
	mean: the rough mean of the read length, default=8000
	return: 'samples' is a list.
			The elements in the list are the length of each reading to be generated
	'''
	samples = st.expon.rvs(213.98910256668592,
		6972.5319847131141, size=size, random_state=seed)
	samples = samples*mean/7106.0
	samples = samples.astype(int)
	samples = np.clip(samples, 1, len(genome))
	return samples

# lambda
# the actual length should be multiplied  by 1000
# two mixture gamma distribution with parameters
# first gamma: alpha: 6.3693711, rate: 0.53834893
# second gamma: alpha: 1.67638771, rate: 0.22871401
def draw_mix_gamma_dis(size, mean, seed):
	'''
	mean: the rough mean of the read length, default=8000
	return: 'samples' is a list.
			The elements in the list are the length of each reading to be generated
	'''
	half = int(size/2.0)
	sample_1 = st.gamma.rvs(6.3693711, 0.53834893, size=half,
		random_state=seed)
	sample = st.gamma.rvs(1.67638771, 0.22871401, size=(size-half),
		random_state=seed)
	sample = np.concatenate((sample, sample_1))
	np.random.seed(seed)
	np.random.shuffle(sample)
	sample = sample*mean/4.39
	sample = sample.astype(int)
	sample = np.clip(sample, 1, len(genome))
	return sample

#----------Part2: Loading the sequence files in txt format or fasta format-------------#

'''
This part is used to load the original sequence.
In order not to generate errors, it is strongly recommended that each input file contains only one original sequence
'''

def load_genome(input_file):
	with open(input_file, 'r') as f:
		text = f.read()
		lines = text.splitlines()
	sequence = filter(lambda x: '>' not in x, lines)
	sequence = map(lambda x: x.strip(), sequence)
	sequence = ''.join(sequence)
	return sequence

def get_start_point(length):
	'''
	return: This function returns the starting position of the reads.
			The random seed is set in the function, and the starting position is randomly selected.
	'''
	return np.random.choice(len(genome)-length+1, 1)[0]

#----------Part3: Generate reads and store in a file in txt or fasta format -------------#

'''
Generate the reads that meet the length distribution requirements.
The reads generated from the same raw sequence are stored in one fasta format file.
'''

def sampling_single_lin(pair):
	'''
	Find the position of the read on the sequence.
	'''
	start_point, read_length = pair
	return genome[start_point: start_point+read_length]

def sampling(read_length, seed):
	'''
	return: 'read_list':Each generated read is saved as an element in this list.
	'''
	np.random.seed(seed)
	read_list = list()
	start_point = map(get_start_point, read_length)
	chunk_pair = zip(start_point, read_length)
	read_list = map(sampling_single_lin, chunk_pair)
	read_list = list(read_list)
	return read_list

def save_file(read_list, output_file):
	'''
	The reads generated from the same raw sequence are stored in one fasta format file.
	'''
	# output_file = output_file + '.fasta'
	output_file = output_file
	with open(output_file, 'w') as f:
		for i in range(len(read_list)):
			f.write('>{}\n'.format(i))
			f.write(read_list[i]+'\n')


def check_coverage(coverage,seq_num_n,len_mean=8000):
	'''
	Compare the specified reading number with the corresponding reading number of the coverage,
	then select the appropriate number of readings.
	'''
	seq_num_c = int(coverage * len(genome) / len_mean)
	if seq_num_c > seq_num_n:
		seq_num = seq_num_c
	else:
		seq_num = seq_num_n
	return seq_num

def check_mean_length(file_path):
	with open(file_path, 'r') as f:
		a = f.read()
		l = a.splitlines()
	l = list(filter(lambda x: '>' not in x, l))
	length = list(map(len, l))
	print('The average length is: ', np.average(length))

############################################################################################################################
# ----------------------------------------------- main -----------------------------------------------#
def initialization_parameters():
	#parser = argparse.ArgumentParser()
	parser = argparse.ArgumentParser(description='Sampling read from the \
		input genome')
	parser.add_argument('-i', action='store', dest='input', type=str, required=True, help='The input file path')
	parser.add_argument('-l', action='store', dest='len_mean',default=8000,
		type=float, help='the rough mean of the read length')
	parser.add_argument('-n', action='store', dest='seq_num', required=True,
		type=int, help='the number of output sequence')
	parser.add_argument('-s', action='store', dest='seed', type=int, default=0,
		help='the random seed, for reproducibility')
	parser.add_argument('-c', action='store', dest='coverage', type=float, default=0,
		help='spacify the simulation coverage, the nubmer of read will be \
		calculated. We use the larger one compared with seq_num.')
	parser.add_argument('-d', action='store', dest='dis', default=3,
		type=int, help='choose from the following distribution: 1: beta_distribution, 2: exponential_distribution, 3: mixed_gamma_dis.The default: 3. ')
	parser.add_argument('-o', action='store', dest='output', type=str, required=True,help='The output file for storing qualified sequences.')
	args = parser.parse_args()
	return args

############################################################################################################################
if __name__ == '__main__':
	time_start = time.time()

	args = initialization_parameters()
	genome = load_genome(args.input)
	seq_num_n = args.seq_num

	seq_num = check_coverage(args.coverage, seq_num_n, args.len_mean)
	#print('The most suitable number is %d' % seq_num + '!')

	dis = args.dis
	if dis == 1:
		#print('You chose the beta_distribution.')
		read_length = draw_beta_dis(seq_num, args.len_mean, args.seed)
		read_list = sampling(read_length, args.seed)
		save_file(read_list, args.output)
	elif dis == 2:
		#print('You chose the exponential_distribution.')
		read_length = draw_expon_dis(seq_num, args.len_mean, args.seed)
		read_list = sampling(read_length, args.seed)
		save_file(read_list, args.output)
	elif dis == 3:
		#print('You chose the mixed_gamma_distribution.')
		read_length = draw_mix_gamma_dis(seq_num, args.len_mean, args.seed)
		read_list = sampling(read_length, args.seed)
		save_file(read_list, args.output)
	else:
		print('Wrong distribution.')

	#print('End!!')

	time_end = time.time()  
	time_sum = time_end - time_start  

	#print('Program running time: %d' % time_sum)