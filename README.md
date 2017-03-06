*TEFLoN*
You must download and install these programs to use TEFLoN.
=======

TEFLoN uses illumina sequence data to discover the positions and identity of TEs present in a sample(s).
Input data can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.

## Dependencies

TEFLoN is only dependent on samtools (www.samtools.sourceforge.net) and the Burrows-Wheeler Aligner (www.bio-bwa.sourceforge.net) if you have an existing TE annotation file in bed format.
If you do not have a TE annotation file in bed format, TEFLoN is also dependent on RepeatMasker (www.repeatmasker.org).

## Usage
###Without reference TE annotation
If you do not have an existing TE annotation, first use teflon_prep_no_anno.py to prepare your reference genome for mapping.

```
usage: python teflon_prep_no_anno.py 
    -wd <full path to the working directory> 
    -r <full path to repeatMasker executable> 
    -g <full path to reference_genome.fa>
    -l <full path to repBase.ref>
    -p <prefix for all newly created files>
    -n <number of CPUs>
```

###With reference TE annotation
Coming soon...



