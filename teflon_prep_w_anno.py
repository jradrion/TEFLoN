import argparse, sys, os
import subprocess as sp
import shlex
import glob
import shutil

teflonBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import reheadRepBaseLib as rh
from teflon_scripts import repeatMaskerOut2Bed as rep2bed
from teflon_scripts import repeatMaskerOut2hierarchy as rep2hier
from teflon_scripts import pseudoGenerate as pg
from teflon_scripts import convertPositions as cp

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-a',dest='anno',help='reference TE annotation')
    parser.add_argument('-t',dest='hier',help='reference TE hierarchy')
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-p',dest='pre',help='prefix for all newly created files')
    args = parser.parse_args()

    #Create directory structure
    dirs_to_create=["prep_TF","prep_MP"]
    for i in range(len(dirs_to_create)):
        D=os.path.join(args.wd,dirs_to_create[i])
        if not os.path.exists(D):
            print "Creating directory:",D
            os.makedirs(D)
    prep_MP_DIR=os.path.join(args.wd,"prep_MP")
    prep_TF_DIR=os.path.join(args.wd,"prep_TF")

    #Copy hierarchy file to prep_TF_DIR
    os.system("cp %s %s" %(args.hier, os.path.join(prep_TF_DIR,args.pre+".hier")))

    #Generate reference in pseudospace
    pseudoRefFILE=os.path.join(prep_MP_DIR,args.pre+".pseudo.fa")
    if not os.path.exists(pseudoRefFILE):
        pickle = pg.pseudo_generate_portal(args.anno,args.genome,prep_MP_DIR,prep_TF_DIR,args.pre)
    else:
        print "Reference in pseudospace already exists:", pseudoRefFILE

    #Convert annotation.bed to pseudospace
    cp.convertPositions_portal(args.anno,pickle,os.path.join(prep_TF_DIR,args.pre+".te.pseudo.bed"))

    #Cat pseudoRef RM.annotatedTE.fa and
    mapRef=os.path.join(prep_MP_DIR,args.pre+".mappingRef.fa")
    if not os.path.exists(mapRef):
        print "Concatonating reference and TE sequences"
        os.system("cat %s %s > %s" %(pseudoRefFILE,os.path.join(prep_MP_DIR,args.pre+".annotatedTE.fa"),mapRef))
    else:
        print "Mapping Reference exists:", mapRef
    print "Reference prep complete."
    print "Map reads to mapping reference:",mapRef

if __name__=="__main__":
    main()
