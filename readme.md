# Oxford Nanopore Sequencing Library Preparation Pipeline

## Overview
***

This project is used for library preparation before Oxford nanopore sequencing. 
It supports barcode library construction. 
At the same time, in addition to the library preparation kit provided by Oxford, 
this project also provides an interface, and users can add custom barcode kits according to their needs.

The steps of library preparation are as follows: 
1. loading the raw files.
2. Fragment the original sequence.
3. Select the required sequencing kit, access primers and sequencing adapters.
4. (Optional) Use DeepSimulator [(https://github.com/liyu95/DeepSimulator)](https://github.com/liyu95/DeepSimulator) 
to simulate the Oxford nanopore sequencing process to obtain sequencing datas.

## Install
### Prerequisites
***
If you want to use [DeepSimulator](https://github.com/liyu95/DeepSimulator) to simulate the process of nanopore sequencing, you
need to choose to install 
Anaconda2 [(https://www.anaconda.com/distribution/)](https://www.anaconda.com/distribution/) or Minoconda2 [(https://conda.io/miniconda.html)](https://conda.io/miniconda.html). 
For example, users may download and install the following Anaconda2 package:

    wget https://repo.anaconda.com/archive/Anaconda2-2018.12-Linux-x86_64.sh
    bash Anaconda2-2018.12-Linux-x86_64.sh


### Packages
***
    You need to install a python interpreter[recommended version: python 3.7.1 or above], 
    then run the following command in a terminal to install the necessary packages to enable the program to run:
        pip install -v numpy==1.20.3
        pip install -v scipy==1.6.2
        pip install -v tqdm==4.55.1
        pip install -v biopython==1.74
        pip install -v pathlib==1.0.1

### Download
***
If you want to use this pipeline, you need to download the package.
You can run the following command to complete the download step: 
        
        Deploy the pipeline package:        
        git clone https://github.com/lykaust15/DeepSimulator.git
        cd ./lib_pre_v.1/     

        Download the DeepSimulator package：
        cd ./tools/
        git clone https://github.com/lykaust15/C.git
        cd ./DeepSimulator/
        ./install.sh

## Run Commands


### Example
***
You can run this program by running the following command in a terminal:  
Here are two simple instructions to use this program: 
    
    python lib_make_v.1.py -i ./example/CP024223.1.fasta -k 8
    python lib_make_v.1.py -i ./example/data_for_96_ETEC -k 3

### Use Barcode Kits
***
You can run the following command to use the barcode-kits function: 

    python lib_make_v.1.py -i ./example/data_for_12_ETEC -k 1 -b ./model/barcode_kits/template2/Native_Barcoding_Expansion_1-12_and_13-24.csv

In this pipeline, it has built-in barcode sequences and the corresponding kits for these sequences commonly used 
by [Oxford Nanopore Sequencing](https://community.nanoporetech.com/docs/sequence/sequencing_software/chemistry-technical-document/v/chtd_500_v1_revae_07jul2016).
The barcode kits supported by the pipeline are as follows:

    EXP-NBD103
    EXP-NBD104
    EXP-NBD114
    SQK-RBK004
    EXP-PBC001
    EXP-PBC096
    SQK-PBK004
    SQK-PCB109
    SQK-RPB004
    SQK-RAB204  
    SQK-16S024

The barcode sequences provided by this pipeline can be found in the following folder:

    ./lib_pre_v.1/model/barcode_kits/template2/
    
This project allows users to provide custom barcodes and flanking sequences on demand. You can run the following command to use this function:

    python lib_make_v.1.py -i ./example/data_for_12_ETEC -f ./model/flanking_template/flanking2.fasta -k 7 -b ./model/barcode_kits/User_template.csv
    
    You will need to provide your flanking sequences and barcode sequences in the following format:
    # -f: the input flanking sequences file, 'FASTA' format. We recommend that users place the flanking sequence files in the folder: ./model/flanking_template/
    we recommend that the use the file name like 'XXXXX.fasta'and the content format as follows: 
    
    >0
    The upstream flanking region for the front part.
    >1
    The upstream flanking region for the rear part.
    >2
    The upstream flanking region for the rear part.
    >3
    The downstream flanking region for the rear part.


    If users only need to add the barcode on one single end, please provide the flanking sequence in the following order:
    
    >0
    The upstream flanking region.
    >1
    The downstream flanking region.
    
    # -k 7: flag to enable custom barcode feature
      -b: the input barcode sequences file, 'CSV' format, We recommend that users place the barcode sequences file in the folder: ./model/barcode_kits/
    we recommend that the use the file name like 'XXXXX_user_barcodes.csv' and the content format as follows: 
    
    |Seq_ID|Reverse Sequence|Forward Sequence|
    |XXXXXX|X..............X|XXXX...........X|
    .....    .......    ......      ...... ...
    |YYYYYY|Y..............Y|YYYY...........Y|
    
    # You can get the CSV file you need by modifying the template file located in this folder: ./model/barcode_kits/User_template.csv
    # You just need to replace the last column in the template with your own barcode sequences.

### Adapter
***
This pipeline also supports user-custom adapter sequences. Users can run the following command to run this program:
    
    -i ./example/data_for_12_ETEC -k 1 -b ./model/barcode_kits/template2/Native_Barcoding_Expansion_1-12_and_13-24.csv -a ./model/adapter_template/flanking.fasta
    
    You will need to provide your own adapter sequences in the following format:
    # -a: the input adapter sequences file, 'FASTA' format. We recommend that users place the adapter sequence files in the folder: ./model/adapter_template/
    we recommend that the use the file name like 'XXXXX.fasta'and the content format as follows: 

    >0
    The top region for the adapter sequence.
    >1 
    The bottom region for the adapter sequence. (Optional)
### Run Deepsimulator
***
You can run this program by running the following command to run the process of simulating nanopore sequencing. 
So that you can get the datas after nanopore sequencing, including electrical signal files (*saved in 'txt' format*) and 
sequencing results (*saved in 'FASTQ' format*):
    
    python ./to_nanopore_seq.py -i ./XXXXXXX/ (input files path) -o ./XXXXXXX/ (output files path)  

The final sequencing results will be saved in this folder _'XXXX_lib_fastq'_ in 
'FASTQ' format. The sequencing results for each raw sequence are stored in a file.
If you want to use more functions of [DeepSimulator](https://github.com/liyu95/DeepSimulator) or change the sequencing 
parameters, you can modify this script: _./ to_nanopore_seq.py_ according to your needs.

### Explanation of the Content in the Output Folder
***
Within the output folder, there are several folders and files. If you run the program with the parameter L turned off, 
then, within the folder _'XXXXXXPre_Flies'_, there are four folders: _'results_loading'_, _'dis_create'_, _'barcode_kits'_, _'lib_results'_. 
Let us explain all of them in chronological order. _'results_loading'_ contains all the sequences that will be used for library preparation. One sequence is stored in each file. 
_'dis_create'_ contains the fragmented results of the raw sequences. All reads from fragmentation of each raw sequence are stored in a file.
_'barcode_kits'_ contains the results of each reading connected to the flanking sequences and barcode sequences. 
_'lib_results'_ contains the final library building preparation (each read in the previous folder accesses a sequencing adapter),
with reads generated from each raw sequence stored in a file. If you run the program with the parameter L turned on, you won't get the _'results_loading'_  folder,
## Parameters
***
For detailed parameter description, please refer to _Parameters.pdf_. It is stored in the folder 'model'. 



