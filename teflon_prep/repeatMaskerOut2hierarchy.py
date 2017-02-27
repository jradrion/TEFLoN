##Converts RepeatMasker out file to a bed file to be used with TEFLoN
import sys, argparse
parser=argparse.ArgumentParser()
parser.add_argument("-b",dest="bed",help="bed file to create")
parser.add_argument("-l",dest="hier",help="output hierarchy file")
args=parser.parse_args()

ID,fam,superfam=[],[],[]
with open(args.bed, "r") as fIN:
    for line in fIN:
        ar = line.split()
        ID.append(ar[3])
        fam.append(ar[3].split("//")[1].split("|")[0])
        superfam.append(ar[3].split("|")[1])

with open(args.hier, "w") as fOUT:
    fOUT.write("%s\t%s\t%s\n" %("id","family","superfamily"))
    for i in range(len(ID)):
        fOUT.write("%s\t%s\t%s\n" %(ID[i],ID[i].split("//")[1].split("|")[0],ID[i].split("|")[1]))


