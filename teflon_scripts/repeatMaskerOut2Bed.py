##Converts RepeatMasker out file to a bed file to be used with TEFLoN
import sys, argparse
parser=argparse.ArgumentParser()
parser.add_argument("-r",dest="rep",help="RepeatMasker out file")
parser.add_argument("-b",dest="bed",help="bed file to create")
args=parser.parse_args()

repMasked=[]
with open(args.rep, "r") as fIN:
    for line in fIN:
        ar=line.split()
        if ar:
            if ar[0] not in "SWscore":
                repMasked.append(ar)

with open(args.bed, "w") as fOUT:
    for line in repMasked:
        if line[8] == "C":
            strand = "-"
        else:
            strand = "+"
        fOUT.write("%s\t%s\t%s\t%s\t.\t%s\n" %(line[4],line[5],line[6],"TEid"+line[14]+"//"+line[9],strand))

