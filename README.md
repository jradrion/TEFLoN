*TEFLoN*
=======

TEFLoN uses illumina sequence data to discover the positions and identity of TEs present in a sample(s).
Input data can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.

## Dependencies

If you have an existing BED formated TE annotation file, TEFLoN is only dependent on samtools (www.samtools.sourceforge.net) and the Burrows-Wheeler Aligner (www.bio-bwa.sourceforge.net) if you have an existing TE annotation file in bed format.
If you do not have a TE annotation file in bed format, TEFLoN is also dependent on RepeatMasker (www.repeatmasker.org).

## Usage
###Without reference TE annotation
Step 1)
If you do not have an existing TE annotation, first use teflon_prep_no_anno.py to prepare your reference genome for mapping.

```
usage: python teflon_prep_no_anno.py 
    -wd <full path to working directory> 
    -r <full path to repeatMasker executable> 
    -g <full path to reference_genome.fa>
    -l <full path to repBase_library.ref>
    -p <prefix for all newly created files>
    -n <number of threads>
```

Step 2)
First, index the prepped reference generated by teflon_prep_no_anno.py
```
bwa index <prefix.mappingRef.fa>
```
Next, map your reads with BWA MEM or a similar mapping package (i.e. one that is able to identify and mark soft-clipped reads).
Note: TEFLoN has only been tested using bwa mem v.0.7.10
```
bwa mem -t <nThreads> -Y <read1.fq> <read2.fq> <alignment.sam>
```
Next, sort and index the alignment
```
samtools view -Sb <alignment.sam> | samtools sort -@ <nThreads> - -o <alignment.sorted.bam>
samtools index <alignment.sorted.bam>
```
Ideally, you should have also QC checked your reads and removed duplicates from the alignment.

Now you are ready to proceed to using TEFLoN

###With reference TE annotation
Coming soon...



