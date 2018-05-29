#########################################################################
# File Name: example.sh                                                 #
# Author: xlzh                                                          #
# mail: xiaolongzhang2015@163.com                                       #
# NOTE: This is an example pipline, which only describe the processing  #
#                                                                       #
# Required:                                                             #
#       1. IDS_TagToHeader.py                                           #
#       2. bwa                                                          #
#       3. samtools                                                     #
#       4. IDS_SSCSMaker.py                                             #
#       5. pilupfilt                                                    #
#       6. Varscan                                                      #
#########################################################################

#!/bin/bash

Data="Sample"
Refer="Refer/Ecoli_ATCC8739"

# move the barcode to header
python IDS_TagToHeader.py \
    --infile1 ${Data}"_R1.fq" \
    --infile2 ${Data}"_R2.fq" \
    --outfile1 ${Data}"_R1.fq.smi" \
    --outfile2 ${Data}"_R2.fq.smi"

# maping
bwa aln ${Refer} ${Data}"_R1.fq.smi" > ${Data}"_R1.aln"
bwa aln ${Refer} ${Data}"_R2.fq.smi" > ${Data}"_R2.aln"
bwa sampe -s ${Refer} ${Data}"_R1.aln" ${Data}"_R2.aln" \
             ${Data}"_R1.fq.smi" ${Data}"_R2.fq.smi" > ${Data}".pe.sam"

# convert the samfile to bamfile
samtools view -Sbu ${Data}".pe.sam" | samtools sort -m 2G -@ 8 -o ${Data}".sort.bam"

# creat single strand consensus sequence (SSCS)
python IDS_SSCSMaker.py \
    -i ${Data}".sort.bam" \
    -t ${Data}".pe.tagcounts" \
    -o ${Data}".sscs.bam" \
    -p 'dpm'

# Creat the mpileup file
samtools mpileup \
    -d 1000000 \
    -Q 20 \
    -o ${Data}".mp"

# filt off the unreliable bases from pileup file
pileupfile ${Data}".mp" ${Data}".filt.mp"

# call snvs
java -jar VarScan.v2.3.7.jar pileup2snp \
     --min-coverage 8 \
     --min-reads2 0 \
     --min-var-freq 0 \
     --p-value 1 > ${Data}".snv"

