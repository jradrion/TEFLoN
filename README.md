*TEFLoN*
=======

TEFLoN uses paired-end illumina sequence data to discover the positions and identity of TEs present in a sample(s).
Input data can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.

## Dependencies

TEFLoN is dependent on samtools (www.samtools.sourceforge.net) and the Burrows-Wheeler Aligner (www.bio-bwa.sourceforge.net) if you have an existing TE annotation file in BED format.
If you do not have a TE annotation file in BED format, TEFLoN is also dependent on RepeatMasker (www.repeatmasker.org).

##Testing
Sample files are provided to ensure TEFLoN and its dependencies are running correctly.

## Usage
###Data prep without reference TE annotation
Step 1) If you do not have an existing TE annotation in BED format, use teflon_prep_no_anno.py to prepare your reference genome for mapping.

```
usage: python /usr/local/teflon_prep_no_anno.py 
    -wd <full path to working directory> 
    -r <full path to repeatMasker executable> 
    -g <full path to reference_genome.fa>
    -l <full path to repBase_library.ref>
    -p <prefix for all newly created files>
    -n <number of threads>
```
Important output files:
* /usr/local/prep_TF/prefix.bed
* /usr/local/prep_TF/prefix.genomeSize.txt
* /usr/local/prep_TF/prefix.hier
* /usr/local/prep_TF/prefix.pseudo2refMap.txt
* /usr/local/prep_TF/prefix.ref2pseudoMap.txt
* /usr/local/prep_MP/prefix.mappingRef.fa

Step 2) Index the prepped reference genome. 
```
bwa index <prefix.mappingRef.fa>
```
Step 3) For each of sample, map reads with BWA mem or a similar mapping software package (i.e. one that is able to identify and mark soft-clipped reads).
Note: TEFLoN has only been tested using bwa mem v.0.7.10
```
bwa mem -t <nThreads> -Y <prefix.mappingRef.fa> <read1.fq> <read2.fq> > alignment.sam
```
Step 4) For each sample, sort and then index the alignment produced by BWA.
```
samtools view -Sb <alignment.sam> | samtools sort -@ <nThreads> - -o alignment.sorted.bam

samtools index <alignment.sorted.bam>
```
Note: Ideally, you should have also QC processed your raw reads and removed duplicates from the alignment.

Step 5) Create a txt file (e.g. samples.txt), where each new line contains both the full path to the indexed and sorted alignment.bam for a sample and a unique prefix/nickname for that sample separated by a tab.
Currently, this step is necessary even if you are only analyzing a single sample.
```
/usr/local/sample1.sorted.bam  s1
/usr/local/sample2.sorted.bam  s2
```

Now you are ready to proceed to using TEFLoN.

###Data prep with reference TE annotation
Coming soon...

###Using TEFLoN
There are four modules to TEFLoN: teflonDiscover, teflonCollapse, teflonCount, and teflonGenotype.
You must run teflonDiscover and teflonCount seperately for each sample. These modules may run independently for each sample (i.e. you can run all samples simultaneosly with if you have enough CPUs).
TeflonCollapse and teflonGenotype only need to run once per analysis.

Step 1) For each sample, run teflonDiscover.
```
usage: python usr/local/teflon.v0.2.py
    -wd <full path to the working directory>
    -ex <full path to the samtools executable>
    -g <full path to genomeSize.txt> #created by teflon_prep
    -s <full path to samples.txt> #user created
    -b <full path to sorted and indexed bam file>
    -a <full path to BED formatted TE annotation file in pseudospace> #created by teflon_prep
    -t <full path to TE hierarchy file> #created by teflon_prep
    -l <level of the hierarchy file to guide initial TE search> #recommended "family" (note: level must appear in the first line of the TE hierarch file)
    -cl <level of the hierarchy to cluster similar TEs> #can be same "level" of hierarchy used in -l or higher (raising level will reduce the number of TE instances found)
    -e <newline separated file of any te families to ignore from analysis> (optional)
    -q <int> #mapped reads with map qualities lower than this number will be discarded
    -sd <int> #use to manually override the insert size sd identified by samtools stat (check this number to ensure it seems more or less correct based on knowledge of sequencing library!)
    -x <int> #number of threads to use
```

Step 2) Run teflonCollapse.

Step 3) For each sample, run teflonCount.

Step 2) Run teflonGenotype.






