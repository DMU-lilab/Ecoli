# Ecoli
Algorithm of Improved Duplex Sequencing

Pipline
====================
### Used to obtain all potential SNVs

    flow:

    IDS_TagToHeader.py
            | (extract barcode to fastq seqname)
            V
         bwa aln
            | (mapping to the reference)
            V
    samtools view
            | (convert to bam format)
            V
    IDS_SSCSMaker.py
            | (creat single strand consensus sequence)
            V
    samtools mpileup
            |
            V
    VarScan.v2.3.7.jar
            | (call snvs with lowest standard)
            V
    All Potential SNVs


SnvFilt
====================
### Get high confidence SNVs with binomial distribution algorithm

    input1: error_rate.txt (12 kinds of error type)

    eg.
    AC      0.0003546
    GT      0.0004029
    AG      0.0005264
    CA      0.0005619
    CG      0.0002339
    GC      0.0002394
    AT      0.0003388
    GA      0.0004577
    CT      0.0004425
    TG      0.0004726
    TC      0.0005198
    TA      0.0003793


    input2: Fulltable.txt (Tab delimited)

    eg.
    Chrom Pos Ref sam1_alt sam1_altnum sam1_cover sam2_alt sam2_altnum sam2_cover ...
    Ecoli 1   G   A        3           235        T        1           233        ...
    Ecoli 2   T   C        1           364        A        2           125        ...
    Ecoli 3   T   G        2           255        C        1           155        ...


Annotation
======================
### Annotate the SNVs

NOTE: The SNV file was obtained by SnvFilt.R

Usage:
    
    python SNPAnno.py <annotation.txt> <SNPfile>


Hotspot
======================
### Get Hotspot region from SNV file

    input: the position of SNV (high confidence, eg. pvalue < 0.01)

    eg.
    Pos     distance
    3471    881
    5169    1698
    6192    1023
    6747    555
    7517    770
    ...     ...

    Note: 1698 = 5169 - 3471; 1023 = 6192 - 5169; ...

Simulate
=======================
### Randomly generate a specified number of mutations

Compile:
    
    gcc -std=c99 -o getsnv

Usage:
    
    ./getsnv <reference.fa> <snv number>


