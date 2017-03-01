*TEFLoN*
=======

TEFLoN uses illumina sequence data to discover the positions and identity of TEs present in a sample/samples.
Input data can either be individually sequenced or pooled data and multiple samples can be analyzed simultaneously to improve sensitivity.

## Dependencies

TEFLoN is only dependent on samtools if you have an existing TE annotation file in bed format.
If you do not have a TE annotation file in bed format, TEFLoN is also dependent on RepeatMasker.

## Usage
If you do not have an existing TE annotation first use teflon_prep_ref.py to prepare your reference genome for mapping

###Without reference TE annotation
First, you will need to use repeatMasker to search for TE sequence in the reference 
```usage: python teflon_prep_ref.py 
    -wd <full path to the working directory> 
    -r <full path to repeatMasker executable> 
    -g <reference_genome.fa>
```

###With reference TE annotation
Coming soon...



