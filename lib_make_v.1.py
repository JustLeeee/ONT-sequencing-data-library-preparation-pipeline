from pack.seq_processing import *
import numpy as np
import pandas as pd
import pathlib
import os
import time
import datetime
import sys
import argparse

############################################################################################################################
# ----------------------------------------------- main function -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='inputfile', type=str, required=True, help='The input file path')
    # parser.add_argument('-L', action='store', dest='search', type=int, default=0, help='Directly search all files in the path and use')
    parser.add_argument('-o', action='store', dest='output_file', type=str, default='0',
                        help='The directory path to save the result files. 0: The default output path. Others: User specified output path')
    parser.add_argument('-l', action='store', dest='len_mean', default=8000,
                        type=float, help='the rough mean of the read length')
    parser.add_argument('-n', action='store', dest='seq_num', default=100,
                        type=int, help='the number of output sequence')
    parser.add_argument('-s', action='store', dest='seed', type=int, default=0,
                        help='the random seed, for reproducibility')
    parser.add_argument('-c', action='store', dest='coverage', type=int, default=0,
                        help='spacify the simulation coverage, the nubmer of read will be calculated. We use the larger one compared with seq_num.')
    parser.add_argument('-d', action='store', dest='dis', default=3,type=int,
                        help='choose from the following distribution: 1: beta_distribution, 2: exponential_distribution, 3: mixed_gamma_dis.The default: 3. ')
    parser.add_argument('-b', action='store', dest='barcode', type=str, default='-1', help='Choose barcode kit, \
                        -1: Native Barcoding Expansion 96 is selected by default')

    parser.add_argument('-e', action='store', dest='ending', type=str, default='N', help='Choose whether to add barcode only in one segment, \
                           N: No. Y: Yes')
    parser.add_argument('-f', action='store', dest='flanking', type=str, default='None',
                        help='The path to the user-defined flanking sequence file.')

    parser.add_argument('-R', action='store', dest='RLB12A',type=int, default=0, help='Choose barcode-RLB12A. 0: No. 1: Yes')
    parser.add_argument('-dir', action='store', dest='direction', type=str, default='F', help='sequence direction, F: Forward. R: Reverse')
    #parser.add_argument('-p', action='store', dest='primer', default='0', type=str, help='Choose primers. 0: Use kit default primers. Others: Use user-defined primers')
    parser.add_argument('-k', action='store', dest='kit', type=int, required=True, help='Choose barcode library kit. \
                        1: Native_Barcoding_Expansion. 2: Rapid_Barcoding_Kit. 3: PCR_Barcoding_Expansion. 4: PCR_Barcoding_Kit. \
                        5: Rapid_PCR_Barcoding_Kit. 6: Barcoding_16S_Kit. 7: User_Custom_Kit. 8: No barcode.' )
    parser.add_argument('-P', action='store', dest='PCR', default=-1,type=int, help='Choose PCR processing. -1: No. Others: Number of PCR')
    parser.add_argument('-A', action='store', dest='Poly_A', type=int, default=-1, help='Choose how to add the A tail: -1: Default---Only one base is added.\
                            The other is the length of adding A tail')
    parser.add_argument('-T', action='store', dest='Poly_T', type=int, default=-1, help='Choose how to add the T tail: -1: Default---Only one base is added.\
                            The other is the length of adding T tail')
    parser.add_argument('-a', action='store', dest='adapter', default='1',type=str,
                        help='Choose an adapter. 1. 1D-adapter. 2. 1D^2-adapter. (3. 2D-adapter). Other: User custom adapter')


    args = parser.parse_args()
    return args

############################################################################################################################
# --------------- main ------------------#
if __name__ == "__main__":
    time_start = time.time()
    print('script name: %s' % sys.argv[0])
    print('script start time: %s\n' % datetime.datetime.now())

    args = initialization_parameters()


    inputfile = args.inputfile
    output_file = args.output_file
    #search = args.search

    len_mean = args.len_mean
    seq_num = args.seq_num
    seed = args.seed
    coverage = args.coverage
    dis = args.dis

    barcode = args.barcode
    RLB12A = args.RLB12A
    direction = args.direction
    #primer = args.primer
    kit = args.kit
    PCR = args.PCR
    ending = args.ending
    flanking = args.flanking

    Poly_A = args.Poly_A
    Poly_T = args.Poly_T
    adapter = args.adapter

    # ---------------------------------------------Step1: Loading the sequences----------------------------------------------#
    file_path = inputfile
    file_path = str(pathlib.Path(file_path))
    # file_path = os.path.splitext(file_path)[0]
    file_name = file_path.split('/')[-1]
    file_name = os.path.splitext(file_name)[0]

    if output_file == '0':
        oo = file_name + 'Pre_Flies'
        outputfile = oo + '/results_loading'
        outputfile = str(pathlib.Path(outputfile))
    else:
        oo = output_file
        outputfile = oo + '/results_loading'
        outputfile = str(pathlib.Path(outputfile))

    f = pathlib.Path(file_path)

    if f.is_file():
        seq_loading_py = "python ./src/seq_loading.py"
        command1_input = ' -i ' + inputfile
        command1_output = ' -o ' + outputfile
        command1_load = seq_loading_py + command1_input + command1_output
        os.makedirs(outputfile)

        print('#---------------------------------Step1------------------------------------#')
        print('Start the loading processing!!')

        os.system(command1_load)
        dirs_path = os.listdir(outputfile)
        aa = len(dirs_path)
        print('Found {} fasta files.'.format(aa))
        print('Finish the loading processing!!!\n')

    if f.is_dir():
        outputfile = inputfile
        dirs_path = os.listdir(outputfile)
        aa = len(dirs_path)
        print('#---------------------------------Step1------------------------------------#')
        print('Start the loading processing!!')
        print('Found {} fasta files.'.format(aa))
        print('Finish the loading processing!!!\n')


    # ---------------------------------------------Step2: Dis_create of the sequences----------------------------------------------#
    dis_create_py = 'python ./src/dis_create.py'
    dirs_path = os.listdir(outputfile)
    outputfile2 = oo + '/dis_create'
    os.makedirs(outputfile2)

    print('#---------------------------------Step2------------------------------------#')
    print('Start the dis_creating processing!!')

    for dir_path in dirs_path:
        num = dir_path.split('_')[-1]
        command2_input = ' -i ' + outputfile + '/' + dir_path
        command2_else= ' -l %f -n %d -s %d -c %f -d %d '%(len_mean,seq_num,seed,coverage,dis)
        command2_output = ' -o ' + outputfile2 + '/dis_create_seq_' + dir_path.split('_')[-1]
        command2_dis = dis_create_py + command2_input + command2_else + command2_output
        os.system(command2_dis)

    print('Finish the dis_creating processing!!!\n')

    # ---------------------------------------------Step3: Choose barcodes and kits----------------------------------------------#
    barcodes_choose_py = 'python ./src/barcodes_choose.py'

    print('#---------------------------------Step3------------------------------------#')
    print('Start choosing barcodes and kits!!\n')


    outputfile3 = oo + '/barcode_kits'
    os.mkdir(outputfile3)

    command3_input = ' -i ' + outputfile2
    #command3_else = ' -b %s -R %d -d %s -p %s -k %d -P %d  -e %s -f %s' %(barcode, RLB12A, direction, primer, kit, PCR, ending, flanking)
    command3_else = ' -b %s -R %d -d %s -k %d -P %d  -e %s -f %s' %(barcode, RLB12A, direction, kit, PCR, ending, flanking)
    command3_output = ' -o ' + outputfile3
    command3_barcodes = barcodes_choose_py + command3_input + command3_else + command3_output
    os.system(command3_barcodes)

    #print(barcode)
    print('\nFinish choosing barcodes and kits!!\n')

    # ---------------------------------------------Step4: Add adapters----------------------------------------------#
    adapter_ligation_py = 'python ./src/adapter_ligation.py'
    dirs_path = os.listdir(outputfile3)

    print('#---------------------------------Step4------------------------------------#')
    print('Start adding adapters!!')


    outputfile4 = oo + '/lib_results'
    os.mkdir(outputfile4)

    for dir_path in dirs_path:
        num = dir_path.split('_')[-1]
        command4_input = ' -i ' + outputfile3 + '/' + dir_path
        command4_else = ' -A %d -T %d -a %s' %(Poly_A, Poly_T, adapter)
        command4_output = ' -o ' + outputfile4 + '/lib_results_seq_' + dir_path.split('_')[-1]
        command4_adapter = adapter_ligation_py + command4_input + command4_else + command4_output
        os.system(command4_adapter)

    print('Finish adding adapters!!\n \n')

    # ---------------------------------------------End----------------------------------------------#
    print('\nFinish library preparation!!\n')




    time_end = time.time()  
    time_sum = time_end - time_start  
    print('Program running time: %d' % time_sum)