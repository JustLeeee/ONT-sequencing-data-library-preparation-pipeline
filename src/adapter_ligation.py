from seq_processing import *
import numpy as np
import os
import time
import argparse
import datetime
import sys

############################################################################################################################
#------------- ----------------------------Part 1: All adapters---------------------------------------------------#
'''
Ligation order: Poly-A tail is connected after the 3' end. 
                Poly-T tail is connected after the 5' end. 
'''

def AddPoly_A_tail(sequence):
    seq_result = sequence + 'A'
    return seq_result

def AddLongPoly_A_tail(sequence,len_tail=0):
    tail = 'A'*len_tail
    seq_result = sequence + tail
    return seq_result

def AddPoly_T_tail(sequence):
    seq_result = 'T' + sequence
    return seq_result

def AddLongPoly_T_tail(sequence,len_tail=0):
    tail = 'T'*len_tail
    seq_result = tail + sequence
    return seq_result

def Addsticky_ends(sequence):
    seq_PA= AddPoly_A_tail(sequence)
    PT_seq_PA= AddPoly_T_tail(seq_PA)
    return PT_seq_PA

#------------- ----------------------------Part 2: All adapters---------------------------------------------------#
'''
The adapters have the following form: 
Adapter sequences for example_seq:  
Adapter Y: top: 5'-GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT -3' 

Here we give a variety of commonly used different adapter sequences officially provided by ONT.
In addition, the ONT official also supports users to provide their own designed adapters according to their needs. 
We also provide such an interface in this code for users to customize the required adapter sequences.
'''

#------------Customize adapters and primers----------------#
def get_adapter(file_name):
    adapter_list = get_seq_list(file_name)
    return adapter_list
'''
def adapter_user_make(sequence1='',sequence2=''):
    adapter_User_top = sequence1
    adapter_User_bottom = sequence2
    return adapter_User_top , adapter_User_bottom'''

def adapter_user(sequence,adapter_User_top='',adapter_User_bottom=''):
    seq_result = adapter_User_top + sequence + adapter_User_bottom
    return seq_result


#------------Ligation Sequencing Kit family----------------#
'''
The adapters in this part are suitable for the 'Ligation Sequencing Kit' family. 
The adapter sequences included in this part are: 
Adapter Mix(AMX), Adapter Mix F(AMX F), Adapter Mix H(AMX H)  
'''

def adapter_mix(sequence):
    # Ligation order: 5'---adapter---3'
    adapter_Y_top = 'GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGC'
    adapter_Y_bottom = 'GCAATACGTAACTGAACGAAGT'
    seq_result = adapter_Y_top + sequence + adapter_Y_bottom
    return seq_result

#------------Rapid Sequencing Kit family----------------#

'''
The adapters in this part are suitable for the 'Rapid Sequencing Kit' family. 
The adapter sequences included in this part are: 
Rapid adapter (RAP)
'''
def rapid_adapter(sequence):
    # Ligation order: 5'---adapter---3'
    adapter_Y_top = 'GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGC'
    adapter_Y_bottom = 'GCAATACGTAACTGAACGAAGT'
    seq_result = adapter_Y_top + sequence + adapter_Y_bottom
    return seq_result

#------------Barcoding and the available kits----------------#
'''
The adapters in this part are suitable for the 'Rapid Sequencing Kit' family. 
This part contains the adapters and primer sequences that need to be added to the kits
'''

def adapter_mix_2(sequence):
    # Ligation order: 5'---adapter---3'
    adapter_Y_top = 'GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGC'
    adapter_Y_bottom = 'GCAATACGTAACTGAACGAAGT'
    seq_result = adapter_Y_top + sequence + adapter_Y_bottom
    return seq_result

#------------RNA and cDNA sequencing and kits----------------#
def Direct_RNA_Sequencing_Kit(sequence,A_tail=10):
    # The Direct RNA Sequencing Kit (SQK-RNA002)
    # RTA
    # Top:
    # 5' - GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC - 3'
    # Bottom:
    # 5' - GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCCTTTTTTTTTT - 3'
    # RMX
    # Top:
    # 5' - TGATGATGAGGGATAGACGATGGTTGTTTCTGTTGGTGCTGATATTGCTTTTTTTTTTTTTATGATGCAAGATACGCAC - 3'
    # Bottom:
    # 5' - GAGGCGAGCGGTCAATTTGCAATATCAGCACCAACAGAAACAACCATCGTCTATCCCTCATCATCAGAACCTACTA - 3'
    seq = AddLongPoly_A_tail(sequence, A_tail) + GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC
    str = 'GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGC'
    adapter_Y_top = str[::-1]
    seq_result = seq + adapter_Y_top
    return seq_result

############################################################################################################################
# ----------------------------------------------- main -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input', type=str, required=True, help='The input file path')
    parser.add_argument('-A', action='store', dest='Poly_A',type=int, default=-1, help='Choose how to add the A tail: -1: Default---Only one base is added.\
                        The other is the length of adding A tail')
    parser.add_argument('-T', action='store', dest='Poly_T',type=int, default=-1, help='Choose how to add the T tail: -1: Default---Only one base is added.\
                        The other is the length of adding T tail')
    parser.add_argument('-a', action='store', dest='adapter', default='1',
                        type=str, help='Choose an adapter. 1. 1D-adapter. 2. 1D^2-adapter. (3. 2D-adapter). Other: User custom adapter')
    parser.add_argument('-o', action='store', dest='output', type=str, required=True,help='The output file for storing qualified sequences.')
    args = parser.parse_args()
    return args

############################################################################################################################
# --------------- main ------------------#
if __name__ == "__main__":
    time_start = time.time()

    '''
    qq = 'NNNATCGATCGNNN'
    #print(adapter_mix(qq))
    #pp = AddLongPoly_A_tail(qq,5)
    pp = AddPoly_A_tail(qq)
    print(pp)
    print(adapter_mix(pp))'''

    #print('script name: %s' % sys.argv[0])
    #print('script start time: %s\n' % datetime.datetime.now())


    args = initialization_parameters()
    seqs_list = get_seq_list(args.input)
    id_list = get_id_list2(args.input)
    results_list = list()
    for seq in seqs_list:
        if args.Poly_A == -1:
            seq = AddPoly_A_tail(seq)
        else:
            seq = AddLongPoly_A_tail(seq,args.Poly_A)

        if args.Poly_T == -1:
            seq = AddPoly_T_tail(seq)
        else:
            seq = AddLongPoly_T_tail(seq,args.Poly_T)

        if args.adapter == '1':
            seq = adapter_mix(seq)
        else:
            l = get_seq_list(args.adapter)
            l = l + ['']
            adapter_seq_top = l[0]
            adapter_seq_bottom = l[1]
            #adapter_seq_top = get_seq_list(args.adapter)[0]
            #adapter_seq_bottom = get_seq_list(args.adapter)[-1]
            #top,bottom = adapter_user_make(adapter_seq_top,adapter_seq_bottom)
            top = adapter_seq_top
            bottom = adapter_seq_bottom
            seq = adapter_user(seq,top,bottom)
        results_list.append(seq)


    a = zip(results_list, id_list)
    #print(results_list)
    in_list = list(a)
    #print(in_list)
    file_name = args.output
    with open(file_name, 'w') as f:
        #print('1')
        #print(len(in_list))
        for line in in_list:
            #print(2)
            f.write(line[1] + '\n')
            #f.write('>' + line[1] + '\n')
            f.write(line[0] + '\n')

    #print('End!!!')

    time_end = time.time()
    time_sum = time_end - time_start
    #print(time_sum)
