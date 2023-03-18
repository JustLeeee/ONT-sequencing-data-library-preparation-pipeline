# -*- coding: utf-8 -*
from seq_processing import *
import numpy as np
import pandas as pd
import os
import time
import datetime
import sys
import argparse

############################################################################################################################
#-------------------------------------------Get barcode kit---------------------------------------------------#
def checkRLB12A(seq_list):
    # check Rapid PCR Barcoding Kit
    seq_list[11] = 'GTTGAGTTACAAAGCACCGATCAG'
    return seq_list

def getbarcodes(file_name,bar=False):
    '''
    Returns a list with all barcodes. 'Native Barcoding Expansion 96 .csv'template is used by default.
    If users need to enter a custom file, the format is as follows:
    No.1(Header): Name1, Name2,......,Barcodes
    No.2: ....................................
    ..........................................
    '''

    if bar:
        data = pd.read_csv(file_name, header=None)
        num_list = list(data.columns)
        seqs_num = num_list[-1]
        barcode_list = list(data[seqs_num])
        barcode_list.pop(0)
        return barcode_list
    else:
        data = pd.read_csv('model/barcode_kits/template2/Native_Barcoding_Expansion_96.csv', header=None)
        num_list = list(data.columns)
        seqs_num = num_list[-1]
        barcode_list = list(data[seqs_num])
        barcode_list.pop(0)
        return barcode_list
'''
def primer_make(file_name):
    l = get_seq_list(file_name)
    l = l + ['']
    F_primer = l[0]
    R_primer = l[1]
    return F_primer , R_primer

def user_primer(F, R, re = False):
    F_primer = F
    R_primer = R
    if re:
        Five_end = R_primer
        Three_end = reverse_selfcomple(F_primer)
        return Five_end, Three_end
    else:
        Five_end = F_primer
        Three_end = reverse_selfcomple(R_primer)
        return Five_end, Three_end
'''
def PCR_process(sequence,n):
    num_amplicon = 2**n
    l = [sequence]*num_amplicon
    return l

#-------------------------------------------Seq_barcode make---------------------------------------------------#
def Native_Barcoding_Expansion(barcode):
    '''
    Use for Native Barcoding Expansion 1-12 (EXP-NBD104) and Native Barcoding Expansion 13-24 (EXP-NBD114):
    5' - AAGGTTAA - barcode - CAGCACCT - 3'
    Native Barcoding Expansion 96(EXP - NBD196):
    Forward primer:
    5' - AAGGTTAA - barcode - CAGCACCT - 3'
    Reverse primer:
    5' - GGTGCTG - barcode - TTAACCTTAGCAAT - 3'
    '''
    #Five_end = 'AAGGTTAA' + barcode + 'CAGCACC'
    Five_end = 'AAGGTTAA' + barcode + 'CAGCACCT'
    Three_end = reverse_selfcomple(Five_end)
    return Five_end, Three_end

def Rapid_Barcoding_Kit(barcode):
    '''
    Use for Rapid Barcoding Kit (SQK-RBK004), Rapid Barcoding Kit 96 (SQK-RBK110.96), Rapid chemistry-based PCR barcoding kits
    5' - GCTTGGGTGTTTAACC - barcode - GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA - 3'
    '''
    Five_end = 'GCTTGGGTGTTTAACC' + barcode + 'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'
    Three_end = reverse_selfcomple(Five_end)
    return Five_end, Three_end

def PCR_Barcoding_Expansion(barcode,re = False):
    '''
    PCR Barcoding Expansion 1-12 (EXP-PBC001) and PCR Barcoding Expansion 1-96 (EXP-PBC096)
    5' - GGTGCTG - barcode - TTAACCT - 3'(full flanking sequence not included)

    The first PCR amplification requires tailed primers to be used which carry these sequences:
    5’ TTTCTGTTGGTGCTGATATTGC-[project-specific forward primer sequence] 3’
    5’ ACTTGCCTGTCGCTCTATCTTC-[project-specific reverse primer sequence] 3’
    '''

    F_primer = 'GGTGCTG' + barcode + 'TTAACCT' + 'TTTCTGTTGGTGCTGATATTGC'
    R_primer = 'GGTGCTG' + barcode + 'TTAACCT' + 'ACTTGCCTGTCGCTCTATCTTC'
    if re:
        Five_end = R_primer
        Three_end = reverse_selfcomple(F_primer)
        return Five_end, Three_end
    else:
        Five_end = F_primer
        Three_end = reverse_selfcomple(R_primer)
        return Five_end, Three_end

def PCR_Barcoding_Kit(barcode,re = False):
    '''
    PCR Barcoding Kit (SQK-PBK004) and the PCR-cDNA Barcoding Kit (SQK-PCB109)
    The top and bottom strand of this primer carry different flanking sequences:
    5' - ATCGCCTACCGTGAC - barcode - ACTTGCCTGTCGCTCTATCTTC - 3'
    5' - ATCGCCTACCGTGAC - barcode - TTTCTGTTGGTGCTGATATTGC - 3'
    The top and bottom sequences are different to avoid 5’ and 3’ end sequences annealing to each other and forming a loop.
    '''
    F_primer = 'ATCGCCTACCGTGAC' + barcode + 'TTTCTGTTGGTGCTGATATTGC'
    R_primer = 'ATCGCCTACCGTGAC' + barcode + 'ACTTGCCTGTCGCTCTATCTTC'

    if re:
        Five_end = R_primer
        Three_end = reverse_selfcomple(F_primer)
        return Five_end, Three_end
    else:
        Five_end = F_primer
        Three_end = reverse_selfcomple(R_primer)
        return Five_end, Three_end

def Rapid_PCR_Barcoding_Kit(barcode):
    '''
    Use for Rapid PCR Barcoding Kit (SQK-RPB004): single rapid attachment primer
    5' - ATCGCCTACCGTGAC - barcode - CGTTTTTCGTGCGCCGCTTC - 3'
    '''
    F_primer = 'ATCGCCTACCGTGAC' + barcode + 'CGTTTTTCGTGCGCCGCTTC'
    R_primer = 'ATCGCCTACCGTGAC' + barcode + 'CGTTTTTCGTGCGCCGCTTC'
    Five_end = F_primer
    Three_end = reverse_selfcomple(R_primer)
    return Five_end, Three_end

def Barcoding_16S_Kit(barcode, base = False, re = False):
    '''
    Use for 16S Barcoding Kit (SQK-RAB204 and SQK-16S024)
    The 3' flanking sequence of the forward primer contains a wobble base (denoted by M; in the primer the base is either an A or a C) in a variable region of the
    16S gene. Default is 'A'
    Forward 16S primer:
    5' - ATCGCCTACCGTGAC - barcode - AGAGTTTGATCMTGGCTCAG - 3'
    Reverse 16S primer:
    5' - ATCGCCTACCGTGAC - barcode - CGGTTACCTTGTTACGACTT - 3'
    '''
    R_primer = 'ATCGCCTACCGTGAC' + barcode + 'CGGTTACCTTGTTACGACTT'
    if base:
        M = 'T'
        F_primer = 'ATCGCCTACCGTGAC' + barcode + 'AGAGTTTGATC' + M + 'TGGCTCAG'
    else:
        M = 'A'
        F_primer = 'ATCGCCTACCGTGAC' + barcode + 'AGAGTTTGATC' + M + 'TGGCTCAG'

    if re:
        Five_end = R_primer
        Three_end = reverse_selfcomple(F_primer)
        return Five_end, Three_end
    else:
        Five_end = F_primer
        Three_end = reverse_selfcomple(R_primer)
        return Five_end, Three_end

def User_Custom_Kit(barcode, flankfile_path, re = False, end = False):
    '''
    Please provide the custom flanking sequence required by the users in the following order:
    The upstream flanking region for the front part.
    The downstream flanking region for the front part.
    The upstream flanking region for the rear part.
    The downstream flanking region for the rear part.
    If users only need to add the barcode on one single end, please provide the flanking sequence in the following order:
    The upstream flanking region.
    The downstream flanking region.
    '''
    if end:
        flank_sequence = get_seq_list(flankfile_path)
        top_flank = flank_sequence[0]
        bottom_flank = flank_sequence[1]

        Five_end = top_flank + barcode + bottom_flank
        Three_end = ''
        return Five_end, Three_end
    else:
        flank_sequence = get_seq_list(flankfile_path)
        F_top_flank = flank_sequence[0]
        F_bottom_flank = flank_sequence[1]
        R_top_flank = flank_sequence[2]
        R_bottom_flank = flank_sequence[3]

        F_primer = F_top_flank + barcode + F_bottom_flank
        R_primer = R_top_flank + barcode + R_bottom_flank

        if re:
            Five_end = R_primer
            Three_end = reverse_selfcomple(F_primer)
            return Five_end, Three_end
        else:
            Five_end = F_primer
            Three_end = reverse_selfcomple(R_primer)
            return Five_end, Three_end

############################################################################################################################
# ----------------------------------------------- main function -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input', type=str, required=True, help='The input file path')
    parser.add_argument('-b', action='store', dest='barcode', type=str, default='-1', help='Choose barcode kit, \
                        -1: Native Barcoding Expansion 96 is selected by default')
    parser.add_argument('-e', action='store', dest='ending', type=str, default='N', help='Choose whether to add barcode only in one segment, \
                        N: No. Y: Yes')
    parser.add_argument('-f', action='store', dest='flanking', type=str, default='None', help='The path to the user-defined flanking sequence file.')
    parser.add_argument('-R', action='store', dest='RLB12A', type=int, default=0, help='Choose barcode-RLB12A. 0: No. 1: Yes')
    parser.add_argument('-d', action='store', dest='direction', type=str, default='F', help='sequence direction, F: Forward. R: Reverse')

    parser.add_argument('-k', action='store', dest='kit', type=int, required=True, help='Choose barcode library kit. \
                        1: Native_Barcoding_Expansion. 2: Rapid_Barcoding_Kit. 3: PCR_Barcoding_Expansion. 4: PCR_Barcoding_Kit. \
                        5: Rapid_PCR_Barcoding_Kit. 6: Barcoding_16S_Kit. 7: User_Custom_Kit. 8: No barcode.')
    parser.add_argument('-P', action='store', dest='PCR', type=int, default=-1, help='Choose PCR processing. -1: No. Others: Number of PCR')
    parser.add_argument('-o', action='store', dest='output', type=str, required=True,help='The output file for storing qualified sequences.')
    args = parser.parse_args()
    return args

############################################################################################################################
# --------------- main ------------------#
if __name__ == "__main__":
    time_start = time.time()

    # ---------------------------------------------Step1: Loading the sequences----------------------------------------------#

    #print('script name: %s' % sys.argv[0])
    #print('script start time: %s\n' % datetime.datetime.now())

    args = initialization_parameters()

    barcode = args.barcode
    #primer = args.primer
    RLB12A = args.RLB12A
    direction = args.direction
    kit = args.kit
    PCR = args.PCR
    ending = args.ending
    flankfile_path = args.flanking

    '''
    if args.primer != '0':
        F_pri, R_pri = primer_make(args.primer)
    else:
        F_pri = ''
        R_pri = ''
    '''
    #raw_path = 'test_for_dis'
    raw_path = args.input  # The path to store all fasta files. 
    dirs_path = os.listdir(raw_path)
    # print(dirs_path)
    l = len(dirs_path) # The number of barcodes to be needed. 

    print('Finish loading the sequences.')
    # ----------------------------------------------Step2: Choose barcodes----------------------------------------------#
    if barcode == '-1':
        bar = False
        barcode_list = getbarcodes(barcode,bar)
        #print(len(barcode_list))
    else:
        bar = True
        barcode_list = getbarcodes(barcode,bar)
        #print(len(barcode_list))

    if RLB12A != 0:
        barcode_list = checkRLB12A(barcode_list)

    print('Finish choose barcodes.')
    # ----------------------------------------------Step3: Choose sequence direction----------------------------------------------#
    if direction == 'F':
        re = False
    elif direction == 'R':
        re = True
    else:
        print('Error direction!!!!!')

    '''
    # noinspection PyUnboundLocalVariable
    user_five,user_three = user_primer(F_pri, R_pri, re)

    print('Finish choose primers')'''
    # ----------------------------------------------Step4: Sequence make----------------------------------------------#

    if ending == 'N':
        end = False
    elif ending == 'Y':
        end = True
    else:
        print('Error ending!!!!!')

    for i in range(0,l):
        #print(dirs_path[i])
        path = raw_path + '/' + dirs_path[i]
        seqs = get_seq_list(path)
        ids = get_id_list(path)
        seq_result = list()
        id_result = list()
        
        if kit == 8:
            for id in ids:
                id = id + '  |  barcode_'  + str(i) + ':  ' + 'None'
                id_result.append(id)
        else:
            for id in ids:
                id = id + '  |  barcode_'  + str(i) + ':  ' + barcode_list[i]
                id_result.append(id)

        for seq in seqs:
            if kit == 1:
                five, three = Native_Barcoding_Expansion(barcode_list[i])
            elif kit == 2:
                five, three = Rapid_Barcoding_Kit(barcode_list[i])
            elif kit == 3:
                five, three = PCR_Barcoding_Expansion(barcode_list[i])
            elif kit == 4:
                five, three = PCR_Barcoding_Kit(barcode_list[i])
            elif kit == 5:
                five, three = Rapid_PCR_Barcoding_Kit(barcode_list[i])
            elif kit == 6:
                five, three = Barcoding_16S_Kit(barcode_list[i])
            elif kit == 7:
                five, three = User_Custom_Kit(barcode_list[i], flankfile_path, re, end)
            elif kit == 8:
                five = ''
                three = ''
            else:
                print('Error barcode!!!!!')
                break

            #seq = user_five + seq + user_three
            seq = five + seq + three
            seq_result.append(seq)

        # ----------------------------------------------Step5: PCR and output----------------------------------------------#
        if PCR == -1:
            a = zip(seq_result, id_result)
            # print(results_list)
            in_list = list(a)
            # print(in_list)
            #file_name = 'result/' + 'T1/' + 'Library_of_' + dirs_path[i]
            file_name = args.output + '/Library_of_' + dirs_path[i]
            #file_name = 'result/' + 'T1' + '/Library_of_' + dirs_path[i]
            with open(file_name, 'w') as f:
                # print('1')
                # print(len(in_list))
                for line in in_list:
                    # print(2)
                    f.write('>' + line[1] + '\n')
                    f.write(line[0] + '\n')
        else:
            k = list()
            n = list()
            for seq in seq_result:
                s = PCR_process(seq,PCR)
                k = k + s

            for id in id_result:
                g = PCR_process(id,PCR)
                n = n + g

            a = zip(k, n)
            # print(results_list)
            in_list = list(a)
            # print(in_list)
            #file_name = 'result/' + 'T1/' + 'Library_of_' + dirs_path[i]
            file_name = args.output + '/Library_of_' + dirs_path[i]
            #file_name = 'result/' + 'T1' + '/Library_of_' + dirs_path[i]
            with open(file_name, 'w') as f:
                # print('1')
                # print(len(in_list))
                h = -1
                for line in in_list:
                    h = h + 1
                    # print(2)
                    f.write('>' + str(h) + '  |Amplicon: seq_' + line[1] + '\n')
                    f.write(line[0] + '\n')

    print('Complete barcodes preparation!!!!')


    time_end = time.time()  
    time_sum = time_end - time_start  
    #print(time_sum)
    #print('Program running time: %d' % time_sum)