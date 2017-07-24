'''generates a pseudo genome by removing all positions specified in the given .bed file
also creates a file of pseudo genome sizes and maps between the pseudo and true reference'''

import os, sys, gzip
import cPickle as pickle

def revComp(sequence):
    '''returns the reverse complement of sequence'''
    comp=""
    for b in sequence:
        if b in "ACGT":
            comp+="TGCA"["ACGT".index(b)]
        else:
            comp+=b
    return comp[::-1]

def fastaformat(seq):
    '''returns a string where every 70 bases is separated by a new line'''
    seq=seq.upper()
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

def removeBedPos(bedFile, fasta, prep_MP_DIR,prep_TF_DIR, pre):
    outFILE=os.path.join(prep_MP_DIR,pre+".pseudo.fa")
    outSEQS=os.path.join(prep_MP_DIR,pre+".annotatedTE.fa")
    bed=[]
    print 'Reading TE annotation:', bedFile
    with open(bedFile, 'r') as fIN:
        for line in fIN:
            arr=line.split()
            chrom=arr[0]
            start=int(arr[1])
            stop=int(arr[2])
            ID=arr[3]
            strand=arr[5]
            bed.append([chrom,start,stop,ID,strand])
    print 'Reading reference:',fasta
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
    extractedSeqs={}
    print "Writing annotated TEs as fasta:",outSEQS
    for ch in genome:
        #print "old len:",len(genome[ch])
        for x in bed:
            if ch == x[0] and x[4]!="-":
                extractedSeqs[x[3]]=""
                for i in xrange(x[1]-1,x[2]):
                    extractedSeqs[x[3]]+=genome[ch][i]
            elif ch== x[0] and x[4]=="-":
                tmp=""
                for i in xrange(x[1]-1,x[2]):
                    tmp+=genome[ch][i]
                extractedSeqs[x[3]]=revComp(tmp)
    for ch in genome:
        for x in bed:
            if ch == x[0]:
                for i in xrange(x[1]-1,x[2]):
                    genome[ch][i]="$"
                #print  extractedSeqs
                #sys.exit()
        pseudoMap[ch]=pseudoMapMake(genome[ch],"$")
        refMap[ch]=refMapMake(genome[ch],"$")
        genome[ch]=remove_values_from_list(genome[ch],"$")
        #print "new len:",len(genome[ch])

    with open(outSEQS, "w") as fOUT:
        for x in extractedSeqs:
            fOUT.write(">%s\n%s\n" %(x,fastaformat(extractedSeqs[x])))

    with open(outFILE, "w") as fOUT1, open(os.path.join(prep_TF_DIR,pre+".genomeSize.txt"), 'w') as fOUT2:
        print "Writing genome size file:",os.path.join(prep_TF_DIR,pre+".genomeSize.txt")
        print "Writing pseudospace to refspace conversion map:",os.path.join(prep_TF_DIR,pre+".pseudo2ref.pickle.gz")
        print "Writing refspace to pseudospace conversion map:",os.path.join(prep_TF_DIR,pre+".ref2pseudo.pickle.gz")
        print "Writing reference in pseudospace:",outFILE
        for x in chroms:
            for ch in genome:
                if ch == x:
                    fOUT1.write(">"+ch+"\n")
                    fOUT1.write(fastaformat("".join(genome[ch]))+"\n")
                    fOUT2.write(ch+"\t1\t"+str(len(genome[ch]))+"\n")

    print "Dumping pickles..."
    pickle.dump(pseudoMap, gzip.open(os.path.join(prep_TF_DIR,pre+".pseudo2ref.pickle.gz"), "wb"))
    print os.path.join(prep_TF_DIR,pre+".pseudo2ref.pickle"), "dumped!"
    pickle.dump(refMap, gzip.open(os.path.join(prep_TF_DIR,pre+".ref2pseudo.pickle.gz"), "wb"))
    print os.path.join(prep_TF_DIR,pre+".ref2pseudo2.pickle"), "dumped!"
    return refMap

def pseudo_generate_portal(bed,genome,prep_MP_DIR,prep_TF_DIR,pre):
    return removeBedPos(bed,genome,prep_MP_DIR,prep_TF_DIR,pre)










