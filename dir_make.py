from pack.seq_processing import *
import os
import time
import datetime
import sys
import random
import argparse

############################################################################################################################
# ----------------------------------------------- main -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', dest='inputfile', type=str, required=True, help='The input file path')
    #parser.add_argument('-seq_len', type=int, required=True,help='Specifies the length of the sequence.')
    parser.add_argument('-s', action='store', dest='seed', type=int, default=0,
                        help='the random seed, for reproducibility')
    parser.add_argument('-output', dest='outputfile', type=str, required=True, help='The output file for storing qualified sequences.')
    args = parser.parse_args()
    return args

##############################################################################################################################
if __name__ == "__main__":
    time_start = time.time()

    #print('script name: %s' % sys.argv[0])
    #print('script start time: %s\n' % datetime.datetime.now())
    args = initialization_parameters()

    inputfile = args.inputfile
    outputfile = args.outputfile
    seed = args.seed

    random.seed(seed)

    file_name = inputfile.split('/')[-1]
    file_name = file_name.split('.')[-1]

    out_form = outputfile + '/' + file_name + '_form.fasta'
    out_re = outputfile + '/' + file_name + '_re.fasta'

    seqlist = get_seq_list(inputfile)
    idlist = get_id_list(inputfile)

    l = len(seqlist)
    n = int(l/2)

    print('There are a total of %d reverse sequences.' %(n))

    aa = zip(seqlist, idlist)
    in_list = list(aa)
    re_list = random.sample(in_list, n)
    #print(a)
    #print(len(a))
    form_list = list(set(in_list) - set(re_list))

    with open(out_form, 'w') as f:
        # print('1')
        # print(len(in_list))
        i = -1
        for line in form_list:
            i = i + 1
            # print(2)
            #f.write(line[1] + '\n')
            f.write('>' + str(i) + '  |  dir: form  >' + line[1] + '\n')
            f.write(line[0] + '\n')

    with open(out_re, 'w') as f:
        # print('1')
        # print(len(in_list))
        i = -1
        for line in re_list:
            i = i + 1
            # print(2)
            #f.write(line[1] + '\n')
            f.write('>' + str(i) + '  |  dir: re  >' + line[1] + '\n')
            f.write(line[0] + '\n')



    time_end = time.time()
    time_sum = time_end - time_start

    print('Program running time: %d' % time_sum)