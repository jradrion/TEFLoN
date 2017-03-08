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


def repeatMask(wd,RM,ref,cpu,repLib):
    try:
        print "Copying %s to %s" %(ref,wd)
        ref_dup=os.path.join(wd,"prep_RM",os.path.basename(ref))
        os.system("cp %s %s" %(ref,ref_dup))
        print "Running RepeatMasker..."
        cmd="perl %s -no_is -nolow -norna -pa %s --lib %s %s" %(RM,cpu,repLib,ref_dup)
        print "cmd:",cmd
        p = sp.Popen(shlex.split(cmd),stdout=sp.PIPE, stderr=sp.PIPE)
        pERR=p.communicate()[1]
        if p.returncode != 0:
            print "Error running RepeatMasker"
            print pERR
            sys.exit(1)
    except OSError:
        print "Cannot run RepeatMasker"
        sys.exit(1)
    print "RepeatMasker complete."

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd',dest='wd',help='full path to working directory')
    parser.add_argument('-r',dest='repMask',help='full path to RepeatMasker executable')
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-l',dest='lib',help='repBase library')
    parser.add_argument('-p',dest='pre',help='prefix for all newly created files')
    parser.add_argument('-n',dest='cpu',help='number of CPUs',type=int,default=1)
    args = parser.parse_args()

    #Create directory structure
    dirs_to_create=["prep_RM","prep_TF","prep_MP"]
    for i in range(len(dirs_to_create)):
        D=os.path.join(args.wd,dirs_to_create[i])
        if not os.path.exists(D):
            print "Creating directory:",D
            os.makedirs(D)
    prep_MP_DIR=os.path.join(args.wd,"prep_MP")
    prep_RM_DIR=os.path.join(args.wd,"prep_RM")
    prep_TF_DIR=os.path.join(args.wd,"prep_TF")

    #Rehead the repBase library
    lib_outFILE=os.path.join(prep_RM_DIR,os.path.basename(args.lib).replace(".ref",".rehead.ref"))
    rh.reheadRepBaseLib_portal(args.lib,lib_outFILE)
    #Mask the reference
    masked_faFILE=os.path.join(prep_RM_DIR,os.path.basename(args.genome)+".masked")
    if not os.path.exists(masked_faFILE):
        repeatMask(args.wd,args.repMask,args.genome,args.cpu,lib_outFILE)

    #Create annotation.bed
    RM_bedFILE=os.path.join(prep_TF_DIR,args.pre+".bed")
    rep2bed.rep2bed_portal(masked_faFILE.replace(".masked",".out"),RM_bedFILE)

    #Create hierarchy.txt
    RM_hierFILE=os.path.join(prep_TF_DIR,args.pre+".hier")
    rep2hier.rep2hier_portal(lib_outFILE,RM_bedFILE,RM_hierFILE)
    #Generate reference in pseudospace
    pseudoRefFILE=os.path.join(prep_MP_DIR,args.pre+".pseudo.fa")
    if not os.path.exists(pseudoRefFILE):
        pg.pseudo_generate_portal(RM_bedFILE,masked_faFILE.replace(".masked",""),prep_MP_DIR,prep_TF_DIR,args.pre)
    else:
        print "Reference in pseudospace already exists:", pseudoRefFILE

    #Cat pseudoRef RM.annotatedTE.fa and
    mapRef=os.path.join(prep_MP_DIR,args.pre+".mappingRef.fa")
    if not os.path.exists(mapRef):
        print "Concatonating reference and TE sequences"
        os.system("cat %s %s %s > %s" %(pseudoRefFILE,os.path.join(prep_MP_DIR,args.pre+".RM.annotatedTE.fa"),lib_outFILE,mapRef))
    else:
        print "Mapping Reference exists:", mapRef
    print "Reference prep complete."
    print "Map reads to mapping reference:",mapRef

if __name__=="__main__":
    main()
