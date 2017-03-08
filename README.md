*TEFLoN*
=======

TEFLoN uses paired-end illumina sequence data to discover the positions and identity of TEs present in a sample(s).
Input data can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.

## Dependencies

TEFLoN is dependent on samtools (www.samtools.sourceforge.net) and the Burrows-Wheeler Aligner (www.bio-bwa.sourceforge.net) if you have an existing TE annotation file in BED format.
If you do not have a TE annotation file in BED format, TEFLoN is also dependent on RepeatMasker (www.repeatmasker.org).

## Usage
###Data prep without reference TE annotation
Step 1) If you do not have an existing TE annotation in BED format, use teflon_prep_no_anno.py to prepare your reference genome for mapping.

```
usage: python teflon_prep_no_anno.py 
    -wd <full path to working directory> 
    -r <full path to repeatMasker executable> 
    -g <full path to reference_genome.fa>
    -l <full path to repBase_library.ref>
    -p <prefix for all newly created files>
    -n <number of threads>
```

Step 2) Index the prepped reference genome. 
```
bwa index <prefix.mappingRef.fa>
```
Step 3) For each of sample, map reads with BWA mem or a similar mapping package (i.e. one that is able to identify and mark soft-clipped reads).
Note: TEFLoN has only been tested using bwa mem v.0.7.10
```
bwa mem -t <nThreads> -Y <read1.fq> <read2.fq> <alignment.sam>
```
Step 4) For each sample, sort and index the alignment produced by BWA.
```
samtools view -Sb <alignment.sam> | samtools sort -@ <nThreads> - -o <alignment.sorted.bam>
samtools index <alignment.sorted.bam>
```
Note: Ideally, you should have also QC processed your raw reads and removed duplicates from the alignment.

Step 5) Create a txt file (e.g. samples.txt), where each new line contains both the full path to the indexed and sorted alignment.bam for a sample and a unique prefix/nickname for that sample separated by a tab.
Currently, this step is necessary even if you are only analyzing a single sample.
```
/usr/local/sample1.bam  s1
/usr/local/sample2.bam  s2
```

Now you are ready to proceed to using TEFLoN.

###Data prep with reference TE annotation
Coming soon...

###Using TEFLoN
There are four modules to TEFLoN: teflonDiscover, teflonCollapse, teflonCount, and teflonGenotype.
You must run teflonDiscover and teflonCount seperately for each sample. These modules may run independently for each sample (i.e. you can run all samples simultaneosly with if you have enough CPUs).
TeflonCollapse and teflonGenotype only need to run once per analysis.

Step 1) For each sample, run teflonDiscover.

Step 2) Run teflonCollapse.

Step 3) For each sample, run teflonCount.

Step 2) Run teflonGenotype.






