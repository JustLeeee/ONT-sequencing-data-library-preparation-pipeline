import numpy as np
import pandas as pd
import shutil
import os
import time
import datetime
import sys
import argparse

############################################################################################################################
# ----------------------------------------------- main function -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_file', type=str, required=True, help='The input file path')
    parser.add_argument('-o', action='store', dest='output_path', type=str, required=True, help='The output file path')
    parser.add_argument('-d', action='store', dest='deepsimulator', type=str, default="./tools/DeepSimulator/deep_simulator.sh",
                        help='The path of deepsimulator.')
    args = parser.parse_args()
    return args

##########################################################################################################################################
# ---------------------------------------------------------- main ------------------------------------------------------#
if __name__ == "__main__":

    time_start = time.time()
    print('script name: %s' % sys.argv[0])
    print('script start time: %s\n' % datetime.datetime.now())

    args = initialization_parameters()


    deepsimulator = args.deepsimulator
    input_file = args.input_file
    output_path = args.output_path
    

    # ------------------------------------ Loading the sequencing----------------------------------------#
    command_else = ' -B 2 -G 1 -n -1 '

    dirs_path = os.listdir(input_file)
    name = output_path.split('/')[-1]
    name = name + '_lib_fastq'
    to_fastq = output_path + '/' + name
    os.mkdir(to_fastq)

    for dir_path in dirs_path:
        fasta_path = input_file + '/' + dir_path
        command_input = ' -i ' + fasta_path
        # deep_output = output_path + '/' + dir_path.split('.')[0]
        deep_output = output_path + '/' + os.path.splitext(dir_path)[0]
        os.mkdir(deep_output)
        command_output = ' -o ' + deep_output
        command_deep = deepsimulator + command_input + command_output + command_else
        os.system(command_deep)
        source = deep_output + '/' + 'pass.fastq'
        # destination = to_fastq + '/' + dir_path.split('.')[0] + '.fastq'
        destination = to_fastq + '/' + os.path.splitext(dir_path)[0] + '.fastq'
        shutil.copy(source, destination)

    time_end = time.time()  
    time_sum = time_end - time_start  
    print('\nSequencing Completed!!!!!')
    print('Program running time: %d' % time_sum)


