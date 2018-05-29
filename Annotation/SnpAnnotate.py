#!/usr/bin/python
import sys

def trans_anno(ann,annlst):
#transform annotation file to annlst
    f=open(ann,'r')
    line = f.readline()
    for eachline in f:
        a=eachline.split('\t')
        a7=a[7].strip()
        annlst.append((a[1],a[2],a[3],a[4],a[5],a[6],a7))
    annlst.remove(annlst[0])
    f.close

def dudevide(annlst,los):
#locus los in annotation
    low=0
    high=len(annlst)-1
    mid=(low+high)//2
    while high-low>1:
        if los > int(annlst[mid][0]):
            low=mid
        else:
            high=mid
        mid=(low+high)//2
    if int(annlst[low][0]) <= los <=int(annlst[low][1]):
        return low
    elif int(annlst[high][0]) <= los <= int(annlst[high][1]):
        return high
    else:
        return low,high

def outputfile(snpfile, annlst):
    f=open(snpfile,'r')
    next(f)
    for eachline in f:
        rawsnp=eachline.strip().split('\t')
        pos=dudevide(annlst,int(rawsnp[1]))
        if type(pos)==tuple:
            a=int(pos[0])
            b=int(pos[1])
            start1=annlst[a][0]
            end1=annlst[a][1]
            start2=annlst[b][0]
            end2=annlst[b][1]
            strand="*"
            productType="*"
            locus="*"
            proteinId="*"
            description="*"
        else:
            start1=annlst[pos][0]
            end1=annlst[pos][1]
            start2="*"
            end2="*"
            locus=annlst[pos][4]
            proteinId=annlst[pos][5]
            description=annlst[pos][6]
            strand=annlst[pos][2]
            productType=annlst[pos][3]
        print "%s\t%s\t%s\t%s\t%s:%s:%s:%s\t%s:%s\t%s:%s\t%s\t%s\t%s\t%s\t%s" % (rawsnp[0],rawsnp[1],rawsnp[2],rawsnp[3],rawsnp[4],rawsnp[5],rawsnp[6],rawsnp[7],start1,end1,start2,end2,strand,productType,locus,proteinId,description)
    f.close

def main ():
    args = sys.argv
    if len(args) < 2:
        print 'Usage: python SNPAnno.py <annotation.txt> <SNPfile>'
        sys.exit(-1)
    annlst = []
    trans_anno(args[1], annlst)
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%("Chrom","position","Ref","Var","Cover:Alt:AltFreq:Pvalue","start1:end1","start2:end2","strand","productType","locus","proteinId","description")
    outputfile(args[2], annlst)

if __name__=='__main__':
    main ()
