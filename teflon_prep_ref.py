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
        cmd="perl "+RM+" -no_is -nolow -norna -pa "+cpu+" --lib "+repLib+" "+ref
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
    parser.add_argument('-r',dest='repMask',help='full path to RepeatMasker')
    parser.add_argument('-g',dest='genome',help='reference genome')
    parser.add_argument('-l',dest='lib',help='RepeatMasker library')
    parser.add_argument('-n',dest='cpu',help='number of CPUs to use')
    args = parser.parse_args()

    rh.reheadRepBaseLib_portal(args.lib)
    repeatMask(args.wd,args.repMask,args.genome,args.cpu,args.lib+".rehead")










if __name__=="__main__":
    main()
