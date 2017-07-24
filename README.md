*TEFLoN*
=======

TEFLoN uses paired-end illumina sequence data to discover the positions and identity of TEs present in one or more datasets.
Input data can either be individually sequenced or sequenced as a pool, and multiple samples can be analyzed simultaneously to improve sensitivity.
*TEFLoN.v0.3 is now ~50X faster!*


## Required packages

* [python](www.python.org) v2.7
* [samtools](www.samtools.sourceforge.net) (tested using v.1.3)
* [Burrows-Wheeler Aligner](www.bio-bwa.sourceforge.net) (tested using v.0.7.10)
* [RepeatMasker](www.repeatmasker.org) (Only required if you do not have a reference TE annotation)

## Testing 
Sample files are provided to ensure TEFLoN and its dependencies are running correctly.

## Usage
### Data prep
Step 1a) If you have a reference TE annotation in BED format, use teflon_prep_w_anno.py to prepare your reference genome for mapping.

NOTE: You must also manually create a file where each entry in the reference TE annotation corresponds at least one label (ideally this label would indicate the family or order for each TE instance, but you can use any label you like.)
TEFLoN will cluster or split reads mapping to the locations specified in the annotation by their respective labels from this user created TE hierarchy file.
The first line of this file must include identifying headers (the header for column of individual TE instances must be "id", the other header labels can be anything and each column must be separated by a tab).
A good example of what this TE hierarchy file should look like is provided in the sample files.

```
usage: python /usr/local/teflon_prep_w_anno.py <required> [optional] 
    -wd [full path to working directory]
    -a <full path to reference TE annotation in BED format> 
    -t <full path to user generated TE hierarchy>
    -g <full path to reference genome in fasta format>
    -p <prefix for all newly created files>
```

Step 1b) If you *do not* have an existing TE annotation in BED format, use teflon_prep_no_anno.py to prepare your reference genome for mapping.

NOTE: In this case you will *not* need to manually create a hierarchy file, as TEFLoN will automatically generate one based on the labels provided in the repBase library. However, you may wish to modify this file to include more information than is automatically generated.

```
usage: python /usr/local/teflon_prep_no_anno.py <required> [optional] 
    -wd [full path to working directory]
    -e <full path to RepeatMasker executable> 
    -g <full path to reference genome in fasta format>
    -r <full path to repBase_library.ref for your organism> #NOTE: a custom TE library may be used, but fasta sequence headers must be formatted as ">FAMILYID\tORDERID\n" 
    -p <prefix for all newly created files>
    -l [minimum length for RepeatMasker predicted TE to be reported in final annotation] 
    -s [RepeatMasker predicted TEs from the same family separated by distances less than the splitDist will be combined into a single annotated TE] #Note: it is recommended to use average sequencing read length 
    -d [only those repeats < x percent diverged from the consensus seq will be included in final annotation]
    -t [number of threads]
```

Step 2) Index the prepped reference genome. 
```
bwa index </usr/local/prefix.prep_MP/prefix.mappingRef.fa>
```
Step 3) For each sample, map all reads using BWA mem with option -Y.
```
bwa mem -t <nThreads> -Y </usr/local/prefix.prep_MP/prefix.mappingRef.fa> <read1.fq> <read2.fq> > alignment.sam
```
Step 4) For each sample, sort and index the alignment produced by BWA.
```
samtools view -Sb <alignment.sam> | samtools sort -@ <nThreads> - -o alignment.sorted.bam

samtools index <alignment.sorted.bam>
```
Note: Ideally, you should have also QC processed your raw reads and removed duplicates from the alignment.

Step 5) Create a .txt file (e.g. samples.txt), where each new line contains both the full path to the indexed and sorted alignment.bam files and a unique prefix/nickname for that sample separated by a tab.
Note that this step is necessary even if you are only analyzing a single sample.

```
Example File:
/usr/local/sample1.sorted.bam  s1
/usr/local/sample2.sorted.bam  s2
```

You are now ready to proceed to using TEFLoN.

### Using TEFLoN
There are four modules in TEFLoN: teflon.v0.3, teflon_collapse, teflon_count, and teflon_genotype.
You must run teflon.v0.3 and teflon_count separately for each sample. These modules run independently for each sample (i.e. you can run all samples simultaneously with enough threads).
Teflon_collapse and teflon_genotype only need to run once per analysis.

Step 1) For each sample, run teflon.v0.3
```
usage: python usr/local/teflon.v0.3.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt> #user created
    -i <unique id for this sample> #must match unique id from samples.txt
    -eb <full path to BWA executable>
    -es <full path to samtools executable>
    -l1 <level of the hierarchy file to guide initial TE search> #it is recommended to use the lowest level in the hierarchy file (i.e. "hier_level_1" for data without user-curated hierarchy)
    -l2 <level of the hierarchy to cluster similar TEs> #same level of the hierarchy used in -l or higher (clustering at higher levels will reduce the number of TE instances found and improve accuracy in determining the TE type)
    -q <map quality threshold> #mapped reads with map qualities lower than this number will be discarded
    -exclude [newline separated file containing the name of any TE families to exclude from analysis] #these names must match names from column -l1 from the hierarchy file
    -sd [insert size standard deviation] #used to manually override the insert size sd identified by samtools stat (check this number in the generated stats.txt file to ensure it seems more or less correct based on knowledge of sequencing library!)
    -cov [coverage override] #used to manually override the estimated coverage if you get the error: "Warning: coverage could not be estimated"
    -t [number of threads]
```

Step 2) Run teflon_collapse.
```
usage: python usr/local/teflon_collapse.py <required> [optional]
    -wd [full path to working directory]
    -s <full path to samples.txt>
    -es <full path to samtools executable>
    -n <TEs must be supported by >= n reads in at least one sample>
    -q <map quality threshold> #mapped reads with map qualities lower than this number will be discarded
    -cov [coverage override] #used to manually override the estimated coverage if you get the error: "Warning: coverage could not be estimated"
    -t [number of threads]
```

Step 3) For each sample, run teflonCount.
```
usage: python /usr/local/teflon_count.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt> #user created
    -i <unique id for this sample> #must match unique id from samples.txt
    -eb <full path to BWA executable>
    -es <full path to samtools executable>
    -l2 <level of the hierarchy to cluster similar TEs> #this should be the same as -l2 from step 1
    -q <map quality threshold>
    -t [number of threads]
```

Step 4) Run teflonGenotype.
NOTE: Currently haploid and diploid data types are under construction!
```
usage: python usr/local/teflon_genotype.py <required> [optional]
    -wd [full path to working directory]
    -d <full path to usr/local/prefix.prep_TF/>
    -s <full path to samples.txt> #user created
    -dt <data type> #must be either haploid, diploid, or pooled #Currently, all types must use pooled for read counts
```

### Output
The output "genotype" file is currently a tab-separated file with the following columns:
```
C1: chromosome
C2: 5' breakpoint estimate ("-" if estimate not available)
C3: 3' breakpoint estimate ("-" if estimate not available)
C4: search level id
C5: cluster level id
C6: strand ("." if strand could not be detected)
C7: annotated reference TE identifier ("-" if novel insertion)
C8: 5' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C9: 3' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C10: read count for "presence reads"
C11: read count for "absence reads"
C12: read count for "other reads"
C13: allele frequency/genotype (frequency if pooled data, genotype for haploid/diploid currently under construction)

```

### Citation
If you use TEFLoN in your work please cite:
```
Adrion, J.R., M.J. Song, D.R. Schrider, M.W. Hahn, and S. Schaack. 2017. Genome-wide estimates of transposable element insertion and deletion rates in *Drosophila melanogaster*. Genome Biology and Evolution. https://doi.org/10.1093/gbe/evx050.
```
