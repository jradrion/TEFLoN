'''convert a bed file with coordinates from a true reference genome to the pseudo genome'''
import argparse, sys, os

parser = argparse.ArgumentParser()
parser.add_argument("-b",dest="bedFILE",help="input bedFILE to be converted")
parser.add_argument("-m",dest="mapFILE",help="input mapFILE")
parser.add_argument("-o",dest="outFILE",help="outFILE")
args = parser.parse_args()


def main():
    print "reading map"
    pseudoMap={}
    with open(args.mapFILE, 'r') as fIN:
        for line in fIN:
            if line.startswith(">"):
                chrom=line.replace("\n","").replace(">","")
            else:
                pseudoMap[chrom]=[int(x) for x in line.split(",")]
    print "converting bed"
    with open(args.bedFILE, 'r') as fIN, open(args.outFILE, 'w') as fOUT:
        for line in fIN:
            ls=line.split()
            chrom,start,stop=ls[0],int(ls[1]),int(ls[2])
            for ch in pseudoMap:
                if ch == chrom:
                    ls[1]=pseudoMap[chrom][start-1]
                    ls[2]=pseudoMap[chrom][stop-1]
                    fOUT.write("\t".join([str(x) for x in ls])+"\n")
main()
