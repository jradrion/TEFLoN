import argparse, sys, os
import subprocess as sp
import shlex
import glob
import shutil

teflonBase = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import pseudoGenerate as pg
from teflon_scripts import ref2pseudoConvert as r2pC

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory',default=0)
    parser.add_argument('-a',dest='anno',help='reference TE annotation')
    parser.add_argument('-t',dest='hier',help='reference TE hierarchy')
    parser.add_argument('-f',dest='fasta',help='canonical TE sequence in fasta format',default=-1)
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-p',dest='pre',help='prefix for all newly created files')
    args = parser.parse_args()

    # set prefix
    prefix=args.pre

    # identify current working directory
    if args.wd == 0:
        cwd=os.getcwd()
    else:
        cwd=os.path.abspath(args.wd)

    #Create directory structure
    dirs_to_create=[prefix+".prep_TF",prefix+".prep_MP"]
    for i in range(len(dirs_to_create)):
        D=os.path.join(cwd,dirs_to_create[i])
        if not os.path.exists(D):
            print "Creating directory:",D
            os.makedirs(D)
    prep_MP_DIR=os.path.join(cwd,prefix+".prep_MP")
    prep_TF_DIR=os.path.join(cwd,prefix+".prep_TF")

    #Copy hierarchy file to prep_TF_DIR
    os.system("cp %s %s" %(args.hier, os.path.join(prep_TF_DIR,prefix+".hier")))

    #Generate reference in pseudospace
    pseudoRefFILE=os.path.join(prep_MP_DIR,prefix+".pseudo.fa")
    pickle = pg.pseudo_generate_portal(args.anno,args.genome,prep_MP_DIR,prep_TF_DIR,prefix)

    #Convert annotation.bed to pseudospace
    r2pC.ref2pseudoConvert_portal(args.anno,pickle,os.path.join(prep_TF_DIR,prefix+".te.pseudo.bed"))

    #Cat pseudoRef RM.annotatedTE.fa and
    mapRef=os.path.join(prep_MP_DIR,prefix+".mappingRef.fa")
    if args.fasta != -1:
        canonicalPATH=args.fasta
        print "Concatonating reference and TE sequences"
        os.system("cat %s %s %s > %s" %(pseudoRefFILE,os.path.join(prep_MP_DIR,prefix+".annotatedTE.fa"),canonicalPATH,mapRef))
    else:
        print "Concatonating reference and TE sequences"
        os.system("cat %s %s > %s" %(pseudoRefFILE,os.path.join(prep_MP_DIR,prefix+".annotatedTE.fa"),mapRef))
    print "Reference prep complete."
    print "Map reads to mapping reference:",mapRef

if __name__=="__main__":
    main()
