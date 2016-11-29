'''generates a pseudo genome by removing all positions specified in the given .bed file
also creates a file of pseudo genome sizes and maps between the pseudo and true reference'''

import argparse, sys, os

parser = argparse.ArgumentParser()
parser.add_argument('-b',dest='bed',help='ref bed')
parser.add_argument('-g',dest='genome',help='reference genome')
parser.add_argument('-o',dest='out',help='ouput file')
args = parser.parse_args()

def fastaformat(seq):
    '''returns a string where every 70 bases is separated by a new line'''
    return '\n'.join(seq[i:i+70] for i in xrange(0,len(seq),70))

def remove_values_from_list(the_list, val):
    return [value for value in the_list if value != val]

def pseudoMapMake(the_list, val):
    return [i+1 for i in xrange(len(the_list)) if the_list[i] != val]

def refMapMake(the_list, val):
    outLine=[]
    count=1
    for x in the_list:
        if x != val:
            outLine.append(count)
            count+=1
        else:
            outLine.append(count)
    return outLine

def removeBedPos(bedFile, fasta, outFILE):
    bed=[]
    print 'reading bed'
    with open(bedFile, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            chrom=arr[0]
            start=int(arr[1])
            stop=int(arr[2])
            bed.append([chrom,start,stop])
    print 'reading genome'
    chroms=[]
    genome={}
    pseudoMap={}
    refMap={}
    with open(fasta, 'r') as fIN:
        chrom=""
        for line in fIN:
            if line.startswith('>'):
                if not chrom:
                    chrom=line.replace(">","").replace("\n","")
                    chroms.append(chrom)
                    seq=[]
                else:
                    print chrom, len(seq), "finished"
                    genome[chrom]=seq
                    seq=[]
                    chrom=line.replace('>','').replace('\n','')
                    chroms.append(chrom)
            else:
                seq.extend([x.upper() for x in line.replace('\n','')])
        print chrom, len(seq), "finished"
        genome[chrom]=seq
    chroms=sorted(chroms)
    print "removing bed positions"
    for ch in genome:
        #print "old len:",len(genome[ch])
        for x in bed:
            if ch == x[0]:
                for i in xrange(x[1]-1,x[2]):
                    genome[ch][i]="$"
        pseudoMap[ch]=pseudoMapMake(genome[ch],"$")
        refMap[ch]=refMapMake(genome[ch],"$")
        genome[ch]=remove_values_from_list(genome[ch],"$")
        #print "new len:",len(genome[ch])

    with open(outFILE, "w") as fOUT1, open(outFILE.replace(".fa",".genomeSize.txt"), 'w') as fOUT2, open(outFILE.replace(".fa",".pseudoMap.txt"),"w") as fOUT3, open(outFILE.replace(".fa",".refMap.txt"),"w") as fOUT4:
        print "writing psuedo genome, genome size file, and conversion map"
        for x in chroms:
            for ch in genome:
                if ch == x:
                    fOUT1.write(">"+ch+"\n")
                    fOUT1.write(fastaformat("".join(genome[ch]))+"\n")
                    fOUT2.write(ch+"\t1\t"+str(len(genome[ch]))+"\n")
                    fOUT3.write(">"+ch+"\n")
                    fOUT3.write(",".join([str(x) for x in pseudoMap[ch]])+"\n")
                    fOUT4.write(">"+ch+"\n")
                    fOUT4.write(",".join([str(x) for x in refMap[ch]])+"\n")



removeBedPos(args.bed, args.genome, args.out)










