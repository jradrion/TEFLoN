*TEFLoN*
=======

TEFLoN uses paired-end illumina sequence data to both discover transposable elements (TEs) and perform TE genotyping.
Input samples can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.


## Required packages

* [python](www.python.org) v2.7
* [samtools](www.samtools.sourceforge.net) (tested using v.1.3)
* [Burrows-Wheeler Aligner](www.bio-bwa.sourceforge.net) (tested using v.0.7.10)
* [RepeatMasker](www.repeatmasker.org) (Only required if you do not have a reference TE annotation)

## Testing 
Test files are provided to ensure TEFLoN and its dependencies are running correctly.

## Usage
### Data prep
Step 1a) If you have a reference TE annotation in BED6 format, use teflon_prep_annotation.py to prepare your reference genome for mapping.

IMPORTANT NOTES:
A unique identifier must be used for each entry in the annotation BED file (column 4). An example annotation file is provided in TEFLoN/test_files.
TEFLoN also required a "TE hierarchy" file, which includes a line corresponding to every TE in the reference annotation and at least one label for that TE (ideally this label would indicate the family/order/class for each TE instance, but you can use any label you like.)
TEFLoN will either group or split TEs by their respective labels provided in this file.
The first line of this tab-separated TE hierarchy file must include identifying headers (the header for column one must be "id", but the other header labels are chosen by the user).
A good example of a properly formatted TE hierarchy file is provided in the test files.

```
usage: python /usr/local/teflon_prep_annotation.py <required> [optional] 
    -wd [full path to working directory]
    -a <full path to reference TE annotation in BED6 format> 
    -t <full path to user generated TE hierarchy>
    -f [canonical TE sequences in fasta format]
    -g <full path to reference genome in fasta format>
    -p <prefix for all newly created files>
```

Step 1b) If you *do not* have an existing TE annotation in BED6 format, use teflon_prep_custom.py to prepare your reference genome for mapping. 

IMPORTANT NOTES:
You will need a library of TE sequences (in fasta format) to map your reads against. Repeat libraries are available via RepBase (http://www.girinst.org/repbase/).
Many RepBase libraries include non-TE sequence in addition to TE sequence. If you include non-TEs in your library, TEFLoN will treat these sequences as if they were TEs, and this will confound your results.
Please take the time to inspect your RepBase library and remove any non-TE sequences.
Please format your TE library such that the description line for each sequence is ">some_unique_te_identifier#te_class\n". 
Most RepBase libraries have order/class information in the description line. You can use this label, provided the description line you submit is formatted as above.
If you do not want to include order/class information, it is still required that your fasta file is formatted as above. Simply replace te_class with any string you like.  
When using teflon_prep_custom.py you will *not* need to manually create a TE hierarchy file, as TEFLoN will automatically generate one based on the labels provided in the TE library (i.e. the te_class label in the description line). 
You may wish to manually edit this hierarchy file to include more information than is automatically generated. It is located at /usr/local/prefix.prep_TF/prefix.hier
An example of a properly formatted TE library file is provided in TEFLoN/test_files.

```
usage: python /usr/local/teflon_prep_no_anno.py <required> [optional] 
    -wd [full path to working directory]
    -e <full path to RepeatMasker executable> 
    -g <full path to reference genome in fasta format>
    -l <full path to TE library for your organism> #NOTE: fasta sequence headers must be formatted as ">some_unique_te_identifier#te_class\n" 
    -p <prefix for all newly created files>
    -c [SW cutoff score for RepeatMasker, default=250] 
    -m [minimum length for RepeatMasker-predicted TE to be reported in the final annotation, default=200] 
    -s [RepeatMasker-predicted TEs that share both the same TE family and strand, separated by distances less than the splitDist will be combined into a single TE instance, default=100] #Note: It is recommended to use your read length 
    -d [only those repeats < x percent diverged from the consensus seq will be included in final annotation, default=20]
    -t [number of threads]
```

Step 2) Index the newly prepared reference genome. 
```
bwa index </usr/local/prefix.prep_MP/prefix.mappingRef.fa>
```
Step 3) For each sample in your analysis, map all reads to mappingRef.fa using BWA mem with option -Y.
```
bwa mem -t <nThreads> -Y </usr/local/prefix.prep_MP/prefix.mappingRef.fa> <read1.fq> <read2.fq> > alignment.sam
```
Step 4) For each sample in your analysis, sort and index the alignment produced by BWA.
```
samtools view -Sb <alignment.sam> | samtools sort -@ <nThreads> - -o alignment.sorted.bam

samtools index <alignment.sorted.bam>
```
Note: Ideally, you should have also QC processed your raw reads and removed duplicates from the alignment.

Step 5) Manually create a tab-separated plain text file (e.g. samples.txt), where each new line contains both the full path to the indexed and sorted alignment.bam file and a unique prefix/nickname for that sample.
This file is required, even if you are only analyzing a single sample.

```
Example File:
/usr/local/sample1.sorted.bam  sample1_name
/usr/local/sample2.sorted.bam  sample2_name
```

You are now ready to proceed to using TEFLoN.

### Using TEFLoN
IMPORTANT NOTES ON WORKFLOW:
There are four main functions in TEFLoN: teflon.v0.4.py, teflon_collapse.py, teflon_count.py, and teflon_genotype.py.
It is required that you run teflon.v0.4.py and teflon_count.py separately for each sample in your analysis. These two modules run independently of one another (i.e. you can run all samples simultaneously).
Teflon_collapse.py and teflon_genotype.py only need to run once per analysis.

Step 1) *For each sample*, run teflon.v0.4.py
```
usage: python usr/local/teflon.v0.4.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt> #NOTE: Manually created by user in previous step
    -i <unique id for this sample> #NOTE: This ID *must match* the unique ID from samples.txt
    -eb <full path to BWA executable>
    -es <full path to samtools executable>
    -l1 <level of the hierarchy file to guide initial TE search> #NOTE: It is recommended that you use the lowest level in the hierarchy file (i.e. "family" for data without a user-curated hierarchy)
    -l2 <level of the hierarchy to group similar TEs> #NOTE: This must be either the same level of the hierarchy used in -l1 or a higher level (clustering at higher levels will reduce the number of TE instances found, but improve accuracy for discriminating TE identity)
    -q <map quality threshold> #NOTE: Mapped reads with map qualities lower than this number will be discarded
    -exclude [newline separated file containing the name of any TE families to exclude from analysis] #NOTE: Use same names as in column one of the hierarchy file
    -sd [insert size standard deviation] #NOTE: Used to manually override the insert size StdDev identified by samtools stat (check this number in the generated stats.txt file to ensure it seems more or less correct based on knowledge of sequencing library!)
    -cov [coverage override] #Note: Used to manually override the coverage estimate if you get the error: "Warning: coverage could not be estimated"
    -t [number of threads]
```

Step 2) *Only once*, run teflon_collapse.py
```
usage: python usr/local/teflon_collapse.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt>
    -es <full path to samtools executable>
    -n1 <TEs must be supported by >= n reads in at least one sample>
    -n2 <TEs must be supported by >= n reads summed across all samples>
    -q <map quality threshold> #NOTE: Mapped reads with map qualities lower than this number will be discarded
    -cov [coverage override] #NOTE: Used to manually override the coverage estimate if you get the error: "Warning: coverage could not be estimated"
    -t [number of threads]
```

Step 3) *For each sample*, run teflonCount.py
```
usage: python /usr/local/teflon_count.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt>
    -i <unique id for this sample>
    -eb <full path to BWA executable>
    -es <full path to samtools executable>
    -l2 <level of the hierarchy to cluster similar TEs> #NOTE: this should be the same as -l2 from teflon.v0.4.py
    -q <map quality threshold>
    -t [number of threads]
```

Step 4) *Only once*, run teflonGenotype.py
```
usage: python usr/local/teflon_genotype.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt>
    -lt [sites genotyped as -9 if adjusted read counts lower than than this threshold, default=1]
    -ht [sites genotyped as -9 if adjusted read counts higher than this threshold, default=mean_coverage + 2*STDEV]
    -dt <data type> #must be either haploid, diploid, or pooled #NOTE: haplid/diploid genotyper under construction, all types must use pooled
```

### Output files
Output "genotype" files are generated for every sample in your analysis, and include the following columns:
These files are located in working_dir/genotypes/
```
C1: chromosome
C2: 5' breakpoint estimate ("-" if estimate not available)
C3: 3' breakpoint estimate ("-" if estimate not available)
C4: search level id (Usually TE family)
C5: cluster level id (Usually TE order or class)
C6: strand ("." if strand could not be detected)
C7: reference TE ID ("-" if novel insertion)
C8: 5' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C9: 3' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C10: read count for "presence reads"
C11: read count for "absence reads"
C12: read count for "ambiguous reads"
C13: genotype for every TE (allele frequency for pooled data, present/absent for halploid, present/absent/heterozygous for diploid) #Note: haploid/diploid caller is under construction, use "pooled" for presence/absence read counts
C14: numbered identifier for each TE in the population

```

### Citation
If you use TEFLoN in your work please cite:
```
Adrion, J.R., M.J. Song, D.R. Schrider, M.W. Hahn, and S. Schaack. 2017. Genome-wide estimates of transposable element insertion and deletion rates in *Drosophila melanogaster*. Genome Biology and Evolution. https://doi.org/10.1093/gbe/evx050.
```
