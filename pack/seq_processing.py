import os
import h5py
import numpy as np
from tqdm import tqdm
import sys
import time
import datetime
import argparse
############################################################################################################################
#---------------Load input sequence---------------#
'''
The original files are in txt format.
seq_list: A list of all sequences
id_list: A list of everything after '>' up to the first '|'. (Note: Do not include '>' and '|')
For example_seq:
   file name: multi.fasta
   id: >0 | 1000000 19463748-20463748 of 50818468
   seq: TGATTCGAAACCCA..............
Then will get:
   id: >0 
   seq: TGATTCGAAACCCA..............
'''

def get_seq_list(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x != '', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    a = list(seq_list)
    return a

def get_id_list(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: x.split('|')[0][1:], lines)
    b = list(id_list)
    return b


def get_id_list2(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    #id_list = map(lambda x: x.split('|')[0][1:], lines)
    #b = list(id_list)
    b = list(lines)
    return b

def groupmake(seq_file):
    seq_list = get_seq_list(seq_file)
    id_list = get_id_list(seq_file)
    #group_num = len(barcode_list)
    a = zip(seq_list, id_list)
    in_list = list(a)
    return in_list


def seq_selfcomple(seq):
    selfComple_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    selfComple_seq = ''
    for item in seq:
        selfComple_seq += selfComple_dict[item]
    return selfComple_seq

def reverse_selfcomple(seq):
    selfComple_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    selfComple_seq = ''
    for item in seq:
        selfComple_seq += selfComple_dict[item]
    lst = list(selfComple_seq)
    lst.reverse()
    reverse_seq = ''.join(lst)
    return reverse_seq

def get_selfseq(file_name):
    seqs = get_seq_list(file_name)
    ids = get_id_list(file_name)
    seqs = map(reverse_selfcomple, seqs)
    return seqs

def get_reseq(file_name):
    seqs = get_seq_list(file_name)
    ids = get_id_list(file_name)
    l = map(seq_selfcomple, seqs)
    a = list(l)
    return a

# ----------------------------------------------- main -----------------------------------------------#
def Initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, action='store', dest='input', required=True, help='The input file path')
    #parser.add_argument('-seq_len', type=int, required=True,help='Specifies the length of the sequence.')
    parser.add_argument('-G', type=int, action='store', dest='get', default=2, help='Get complementary and reverse complementary sequences. \
                        1: Get complementary;2: Get complementary and reverse complementary sequences')
    parser.add_argument('-p', type=int, action='store', dest='make_flie', default=0, help='The output file for storing qualified sequences.\
                        0: Do not make a fasta file. 1: Make a fasta file.')
    parser.add_argument('-output', type=str, action='store', dest='output', help='The output file for storing qualified sequences.')
    args = parser.parse_args()
    return args
############################################################################################################################
# --------------- main ------------------#
if __name__ == "__main__":
    time_start = time.time()

    args = Initialization_parameters()

    seqs = get_seq_list(args.input)
    ids = get_id_list(args.input)
    result = list()

    if args.get == 1:
        for seq in seqs:
            s = seq_selfcomple(seq)
            result.append(s)
    elif args.get == 2:
        for seq in seqs:
            s = reverse_selfcomple(seq)
            result.append(s)
    print(result)

    if args.make_file !=0:
        a = zip(result, ids)
        in_list = list(a)
        file_name = args.output
        with open(file_name, 'w') as f:
            for line in in_list:
                f.write('>' + line[1] + '\n')
                f.write(line[0] + '\n')


    print('Finish!!!')


    time_end = time.time()  
    time_sum = time_end - time_start  
    print('Program running time: %d' % time_sum)