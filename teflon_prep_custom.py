import argparse, sys, os
import subprocess as sp
import shlex

teflonBase = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.insert(1, teflonBase)

from teflon_scripts import repeatMaskerOut2bed as rep2bed
from teflon_scripts import repeatMaskerOut2hierarchy as rep2hier
from teflon_scripts import pseudoGenerate as pg
from teflon_scripts import ref2pseudoConvert as r2pC

def repeatMask(wd,prefix,RM,ref,cpu,repLib,div,cutoff):
    try:
        print "Copying %s to %s" %(ref,wd)
        ref_dup=os.path.join(wd,prefix+".prep_RM",os.path.basename(ref))
        os.system("cp %s %s" %(ref,ref_dup))
        print "Running RepeatMasker..."
        cmd="perl %s -a -div %s -no_is -nolow -norna -cutoff %s -pa %s --lib %s %s" %(RM,div,cutoff,cpu,repLib,ref_dup)
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
    parser.add_argument('-wd',dest='wd',help='full path to working directory', default=0)
    parser.add_argument('-e',dest='repMask',help='full path to RepeatMasker executable')
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-l',dest='lib',help='custom library')
    parser.add_argument('-p',dest='pre',help='prefix for all newly created files')
    parser.add_argument('-c',dest='cutoff',help='Sets cutoff score for masking repeats',type=int,default=225)
    parser.add_argument('-m',dest='minLen',help='minimum length for RM predicted TE to be reported in annotation',type=int,default=200)
    parser.add_argument('-s',dest='splitDist',help='RM predicted TEs of the same family separated by distances less than the splitDist will be combined into a single annotated TE',type=int, default=100)
    parser.add_argument('-d',dest='div',help='masks only those repeats < x percent diverged from consensus seq',type=int, default=20)
    parser.add_argument('-t',dest='cpu',help='number of threads',type=int,default=1)
    args = parser.parse_args()

    # set prefix
    prefix=args.pre

    # identify current working directory
    if args.wd == 0:
        cwd=os.getcwd()
    else:
        cwd=os.path.abspath(args.wd)

    #Create directory structure
    dirs_to_create=[prefix+".prep_RM",prefix+".prep_TF",prefix+".prep_MP"]
    for i in range(len(dirs_to_create)):
        D=os.path.join(cwd,dirs_to_create[i])
        if not os.path.exists(D):
            print "Creating directory:",D
            os.makedirs(D)
    prep_MP_DIR=os.path.join(cwd,prefix+".prep_MP")
    prep_RM_DIR=os.path.join(cwd,prefix+".prep_RM")
    prep_TF_DIR=os.path.join(cwd,prefix+".prep_TF")

    #Mask the reference
    masked_faFILE=os.path.join(prep_RM_DIR,os.path.basename(args.genome)+".masked")
    repeatMask(cwd,prefix,args.repMask,args.genome,args.cpu,args.lib,args.div,args.cutoff)

    #Create annotation.bed
    RM_bedFILE=os.path.join(prep_RM_DIR,prefix+".bed")
    rep2bed.rep2bed_portal(masked_faFILE.replace(".masked",".out"),RM_bedFILE,args.splitDist,args.minLen)

    #Create hierarchy.txt
    RM_hierFILE=os.path.join(prep_TF_DIR,prefix+".hier")
    rep2hier.rep2hier_portal(args.lib,RM_bedFILE,RM_hierFILE)
    #Generate reference in pseudospace
    pseudoRefFILE=os.path.join(prep_MP_DIR,prefix+".pseudo.fa")
    pickle = pg.pseudo_generate_portal(RM_bedFILE,masked_faFILE.replace(".masked",""),prep_MP_DIR,prep_TF_DIR,prefix)

    #Convert annotation.bed to pseudospace
    r2pC.ref2pseudoConvert_portal(RM_bedFILE,pickle,os.path.join(prep_TF_DIR,prefix+".te.pseudo.bed"))
    #Cat pseudoRef RM.annotatedTE.fa and
    mapRef=os.path.join(prep_MP_DIR,prefix+".mappingRef.fa")
    print "Concatonating reference and TE sequences"
    os.system("cat %s %s %s > %s" %(pseudoRefFILE,os.path.join(prep_MP_DIR,prefix+".annotatedTE.fa"),args.lib,mapRef))
    print "Reference prep complete."
    print "Map reads to mapping reference:",mapRef

if __name__=="__main__":
    main()
