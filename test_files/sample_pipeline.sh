#install dependencies using the following:
#conda create --name teflon_env --yes -c bioconda bwa samtools=1.3 python=2.7
#conda activate teflon_env

###Define variables
WD="./sample_output/"
PPN="4"
PREFIX="TEST"
HIERARCHY="./sample_reference_hierarchy.txt"
ANNOTATION="./sample_reference_te_annotation.bed"
GENOME="./sample_reference_genome.fasta"
READS1="./sample_reads_1.fq"
READS2="./sample_reads_2.fq"
SAMPLES="./sample_names.txt"


##Prep reference and map reads for TEFLoN
python ./../teflon_prep_annotation.py \
    -wd ${WD}reference/ \
    -a $ANNOTATION \
    -t $HIERARCHY \
    -g $GENOME \
    -p $PREFIX

bwa index ${WD}reference/${PREFIX}.prep_MP/${PREFIX}.mappingRef.fa

bwa mem -t $PPN -Y \
    ${WD}reference/${PREFIX}.prep_MP/${PREFIX}.mappingRef.fa \
    $READS1 \
    $READS2 \
    > ${WD}reference/${PREFIX}.sam

samtools view \
    -Sb ${WD}reference/${PREFIX}.sam \
    | samtools sort \
    -@ $PPN - \
    -o ${WD}reference/${PREFIX}.sorted.bam

samtools index \
    ${WD}reference/${PREFIX}.sorted.bam
rm ${WD}reference/${PREFIX}.sam


##Run TEFLoN
##teflon discover
python ./../teflon.v0.4.py \
    -wd ${WD} \
    -d ${WD}reference/${PREFIX}.prep_TF/ \
    -s ${SAMPLES} \
    -i sample1 \
    -l1 family \
    -l2 order \
    -q 20 \
    -t $PPN

##teflon collapse
python ./../teflon_collapse.py \
    -wd ${WD} \
    -d ${WD}reference/${PREFIX}.prep_TF/ \
    -s ${SAMPLES} \
    -n1 1 \
    -n2 1 \
    -q 20 \
    -t $PPN

##teflon count
python ./../teflon_count.py \
    -wd ${WD} \
    -d ${WD}reference/${PREFIX}.prep_TF/ \
    -s ${SAMPLES} \
    -i sample1 \
    -l2 order \
    -q 20 \
    -t $PPN

##teflon genotype
python ./../teflon_genotype.py \
    -wd ${WD} \
    -d ${WD}reference/${PREFIX}.prep_TF/ \
    -s ${SAMPLES} \
    -lt 1 \
    -ht 100 \
    -dt pooled
