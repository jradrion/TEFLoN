import argparse, sys, os
import multiprocessing as mp
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


def repeatMask(wd,RM,ref,cpu,repLib):
    try:
        print "Copying %s to %s" %(ref,wd)
        ref_dup=os.path.join(wd,"repeatMasker_files",os.path.basename(ref))
        os.system("cp %s %s" %(ref,ref_dup))
        print "Running RepeatMasker..."
        cmd="perl %s -no_is -nolow -norna -pa %s --lib %s %s" %(RM,cpu,repLib,ref_dup)
        print "cmd:",cmd
        p = sp.Popen(shlex.split(cmd),stdout=sp.PIPE, stderr=sp.PIPE)
        pOUT=p.communicate()[0]
        pERR=p.communicate()[1]
        if p.returncode != 0:
            print "Error running RepeatMasker"
            sys.exit(1)
    except OSError:
        print "Cannot run RepeatMasker"
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-r',dest='repMask',help='full path to RepeatMasker executable')
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-l',dest='lib',help='repBase library')
    parser.add_argument('-p',dest='pre',help='prefix for all newly created files')
    parser.add_argument('-n',dest='cpu',help='number of CPUs',type=int,default=1)
    args = parser.parse_args()

    dirs_to_create=["repeatMasker_files","teflon_input_files","mapping_files"]
    for i in range(len(dirs_to_create)):
        D=os.path.join(args.wd,dirs_to_create[i])
        if not os.path.exists(D):
            print "Creating directory:",D
            os.makedirs(D)

    #Rehead the repBase library
    lib_outFILE=os.path.join(args.wd,"repeatMasker_files",os.path.basename(args.lib).replace(".ref",".rehead.ref"))
    rh.reheadRepBaseLib_portal(args.lib,lib_outFILE)

    #Mask the reference
    ref_dup=os.path.join(args.wd,"repeatMasker_files",os.path.basename(args.genome))
    masked_dup=ref_dup+".masked"
    RM_outFILE=ref_dup+".out"
    if not os.path.exists(RM_outFILE):
        repeatMask(args.wd,args.repMask,args.genome,args.cpu,lib_outFILE)

    #Create annotation.bed
    inputDIR=os.path.join(args.wd,"teflon_input_files")
    mappingDIR=os.path.join(args.wd,"repeatMasker_files")
    RM_bedFILE=os.path.join(inputDIR,args.pre+".bed")
    rep2bed.rep2bed_portal(RM_outFILE,RM_bedFILE)

    #Create hierarchy.txt
    RM_hierFILE=os.path.join(inputDIR,args.pre+".hier")
    rep2hier.rep2hier_portal(RM_bedFILE,RM_hierFILE)

    #Generate pseudo genome
    pg.pseudo_generate_portal(RM_bedFILE,masked_dup,mappingDIR,inputDIR,args.pre)












if __name__=="__main__":
    main()
