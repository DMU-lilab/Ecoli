#!/usr/bin/env python
'''
SSCSMarker
Inputs:
    A position-sorted paired-end BAM file containing reads with a duplex tag in the header.
Outputs:
    1. output bam file(After filter the duplicate reads)
    2. SSCS.tagcount file
    3. SSCS.tagstate file

    Special arguments description:
    [-p] [-read_type] : A string specifying which types of read to consider.
                        Read types as follows [default: 'dpm']: 
                        n: Neither read 1 or read 2 mapped. 
                        m: Either read 1 or read 2 mapped, but not both. 
                        p: Both read 1 and read 2 mapped, not a propper pair. 
                        d: Both read 1 and read 2 mapped, propper pair. 
                        s: Single ended reads.
    [-r] [-rep_filt] :  Remove tags with homomeric runs of nucleotides of length x. [default: 9]
'''

import sys
import re
import os.path
import pysam
from collections import defaultdict
from optparse import OptionParser

#global regular to match the tag (eg: "|CCAATGCCGCTGATGCGCAAACTA/") in seqname
TagMatch = re.compile('\|[ATCG_]+/')


def TagGet (read, repcount):
    try:
        if read.is_read1:
            readmark = ":1"
        elif read.is_read2:
            readmark = ":2"
        else:
            readmark = ":se"
        tag = TagMatch.search(read.qname).group()[1:-1] + readmark
        return tag
    except:
        return 'A' * repcount


def GoodFlagDefine (flagtype):
    goodFlag = []
    if 'd' in flagtype:
        goodFlag.extend((99, 83, 163, 147))
    if 'm' in flagtype:
        goodFlag.extend((181, 117, 137, 133, 73, 89, 69, 153))
    if 'p' in flagtype:
        goodFlag.extend((97, 81, 161, 145, 129, 65, 177, 113))
    if 'n' in flagtype:
        goodFlag.extend((141, 77, 4))
    if 's' in flagtype:
        goodFlag.extend((0, 16))
    return goodFlag


def OverlapCheck (read):
    overlap, readlen = False, len(read.seq)
    if read.pos < read.mpos and read.mpos < read.pos + readlen and int(read.flag) in (83, 99, 147, 163):
        overlap = True
    elif read.pos > read.mpos and read.pos < read.mpos + readlen and int(read.flag) in (83, 99, 147, 163):
        overlap = True
    elif read.pos == read.mpos and int(read.flag) in (83, 99, 147, 163):
        overlap = True
    return overlap


def TagFlagLenCheck (read, flagtype, repfilt, tag):
    matchlen ,readlen = 0, len(read.seq)
    checkstatus = True
    goodFlag = GoodFlagDefine(flagtype)
    badTag = ('A'*repfilt, 'T'*repfilt, 'G'*repfilt, 'C'*repfilt)
    for cig in read.cigar:
        if cig[0] == 0:
            matchlen += cig[1]
    if matchlen < readlen * 0.95:
        checkstatus = False
    if (read.flag not in goodFlag) or (tag in badTag):
        checkstatus = False
    return checkstatus


def CreatSSCSRead (readlist):
    consensus = ''
    nucIdentity = {'A':0,'T':0,'G':0,'C':0,'N':0}
    totalcol, readlen = 0, len(readlist[0].seq)
    if len(readlist) == 2:
        for i in xrange(1,readlen):
            try:
                if readlist[0].seq[i] != readlist[1].seq[i]:
                    if readlist[0].qual[i] > readlist[1].qual[i]:
                        return readlist[0]
                    else:
                        return readlist[1]
            except:
                return readlist[0]
    else:
        for i in xrange(readlen):
            for j in xrange(len(readlist)):
                try:
                    nucIdentity[readlist[j].seq[i]] += 1
                    totalcol += 1
                except:
                    break
            maxBase = max(nucIdentity, key=nucIdentity.get)
            try:
                if float(nucIdentity[maxBase]) / float(totalcol) > 0.6:
                    consensus += maxBase
                nucIdentity = {'A':0,'T':0,'G':0,'C':0,'N':0}
                totalcol = 0
            except:
                break 
        for i in xrange(len(readlist)):
            if not re.search(consensus,readlist[i].seq):
                return readlist[i]
            elif not re.search(consensus[1:-1],readlist[i].seq):
                return readlist[i]
        return False
                
        
def ConsensusMaker (readDict, readcount):
    consensusDict = {'singleread': [], 'sscsread':[]}
    for tagkey in readDict.keys():
        if len(readDict[tagkey]) == 1:
            consensusDict['singleread'].append(readDict[tagkey][0])
            readcount['singlefamilynum'] += 1
        elif len(readDict[tagkey]) >= 2:
            sscsread = CreatSSCSRead(readDict[tagkey])
            if sscsread:
                consensusDict['sscsread'].append(sscsread)
                readcount['sscsnum'] += 1
    return consensusDict


def WriteOut (consensusDict, outbam):
    for sscsread in consensusDict['sscsread']:
        outbam.write(sscsread)

    for singleread in consensusDict['singleread']:
        outbam.write(singleread)


def TagInfoWrite (tagDict, tagfile):
    with open (tagfile, "w") as tagcountfile:
        for tagkey in sorted(tagDict.keys()):
            tagcountfile.write('%s\t%d\n' %(tagkey, tagDict[tagkey]))

    totalreads = 0
    familycount = defaultdict ( lambda: 0 )
    with open (tagfile.replace('.tagcounts', '.tagstats'), 'w') as tagstatfile:
        for tagvalue in tagDict.values():
            familycount[tagvalue] += 1
        for size in familycount.keys():
            familycount[size] *= int(size)
            totalreads += int(familycount[size])
        for size in sorted(familycount.keys()):
            tagstatfile.write("%s\t%s\n" %(size, float(familycount[size])/float(totalreads)))


def main():
    usage = 'Usage: %prog -i <input.bam> -o <out.bam> -l <realseq length>'
    parser = OptionParser(usage=usage)
    parser.add_option('-i', "--infile", action="store", dest="infile", 
            help="[required] input BAM file")
    parser.add_option('-o', "--outfile",  action="store", dest="outfile", 
            help="[required] output BAM file")
    parser.add_option('-t', "--tagfile",  action="store",  dest="tagfile", 
            help="output tagcounts file",  default='SSCS.tagcounts')
    parser.add_option('-r', "--rep_filt", action="store",  type=int, dest='rep_filt', 
            help="Remove tags with sam nucleotides of length x. [6]", default=6 )
    parser.add_option('-p', '--read_type', type=str, action="store", dest='read_type', default="dpm", 
            help="Read types: n, m, p, d, s. default = ['dpm']")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit(0)

    InBam = pysam.Samfile(options.infile, "rb")
    OutBam = pysam.Samfile(options.outfile,"wb", template = InBam)
    BamEntry = InBam.fetch(until_eof = True)

    previousread, currentread = BamEntry.next(), ''
    readcount = {'allreadnum':1, 'singlefamilynum':0, 'sscsnum':0}
    readDict, consensusDict, tagDict = {}, {}, defaultdict( lambda: 0 )
    ''' DataStruct (examples as fellows)
        readDict = {tag1:[read1,read2,...], tag2:[read1,read2,...], ...}
        consensusDict = {singleread:[read1,read2,...], sscsread:[read1,read2,...]} (Only have two members)
        tagDict = {tag1:15, tag2: 187, tag3:18, ...}
    '''
    readDict[TagGet(previousread, options.rep_filt)] = [previousread]
    filedone = False

    while (filedone == False):
        try:
            currentread = BamEntry.next()
            readcount['allreadnum'] += 1
        except:
            filedone = True
            consensusDict = ConsensusMaker(readDict, readcount)
            WriteOut (consensusDict, OutBam)
            print 'I am coming ...'
#            TagInfoWrite (tagDict, options.tagfile)
            continue

        if (currentread.pos == previousread.pos):
            if readcount['allreadnum'] % 100000 == 0:
                sys.stdout.write('\r[*]  Reads processed: %d' %(readcount['allreadnum']))
                sys.stdout.flush()
            tag = TagGet(currentread, options.rep_filt)
#            tagDict[tag] += 1
            if (TagFlagLenCheck(currentread, options.read_type, options.rep_filt, tag)\
                and not OverlapCheck(currentread)):
                if tag not in readDict.keys():
                    readDict[tag] = [currentread]
                else:
                    readDict[tag].append(currentread)
            previousread = currentread
        else:
            consensusDict = ConsensusMaker(readDict, readcount)
            WriteOut (consensusDict, OutBam)
            readDict = {}

            tag = TagGet(currentread, options.rep_filt)
#            tagDict[tag] += 1
            if (TagFlagLenCheck(currentread, options.read_type, options.rep_filt, tag) \
                and not OverlapCheck(currentread)):
                readDict[tag] = [currentread]
            previousread = currentread


    # Close the files opened
    InBam.close()
    OutBam.close()

    # Write summary statistics
    sys.stdout.write("\n")
    sys.stdout.write("Summary Statistics: \n")
    sys.stdout.write("Total reads processed : %d\n" %readcount['allreadnum'])
    sys.stdout.write("Reads after filtered : %d\n" %(readcount['singlefamilynum'] + readcount['sscsnum']))
    sys.stdout.write("\t(1)Reads of siglefamily : %d\n" %readcount['singlefamilynum'])
    sys.stdout.write("\t(2)Reads os SSCS : %d\n" %readcount['sscsnum'])


if __name__ == '__main__':
    main()

