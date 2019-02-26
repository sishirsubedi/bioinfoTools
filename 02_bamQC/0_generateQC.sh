#!/usr/bin/env bash
export SHELL=/usr/bin/bash


ENV="test"
REF_1="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
#################################################
# Parsing arguments
#################################################

if [ "$#" -eq 0 ]; then
echo "Usage: runAlignment.sh"
echo "-s Sample Name"
echo "-f Fastq Dir"
echo "-o Output Dir"
exit
fi

while getopts :s:f:o: option; do
	case "$option" in
    s) SAMPLE="$OPTARG" ;;
    f) FASTQ_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    :) echo "Option -$OPTARG requires an argument." ;;
	  \?) echo "Invalid option: -$OPTARG" ;;
	esac
done

OUTPUT_DIR_SAMPLE ="${OUTPUT_DIR}${SAMPLE}/"

if [ ! -d $OUTPUT_DIR_SAMPLE ]
then
  mkdir $OUTPUT_DIR_SAMPLE
fi
chmod 777 $OUTPUT_DIR_SAMPLE

OUTPUT_DIR_SAMPLE_ALIGNMENT ="${OUTPUT_DIR}${SAMPLE}Alignment/"

if [ ! -d $OUTPUT_DIR_SAMPLE_ALIGNMENT ]
then
  mkdir $OUTPUT_DIR_SAMPLE_ALIGNMENT
fi
chmod 777 $OUTPUT_DIR_SAMPLE_ALIGNMENT

echo "Hello" > ${OUTPUT_DIR_SAMPLE_ALIGNMENT}$sample.trimmomatic.summary.txt


# echo "Running Trimmomatic: Removing sequences < Q20 sample- $SAMPLE"
# trimmomatic="sudo java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 15 \
#               $FASTQ_DIR/"$SAMPLE"_R1_001.fastq.gz \
#               $FASTQ_DIR/"$SAMPLE"_R2_001.fastq.gz \
#               $FASTQ_DIR/"$SAMPLE"_filt_paired_q20_R1_001.fastq.gz \
#               $FASTQ_DIR/"$SAMPLE"_filt_unpaired_q20_R1_001.fastq.gz \
#               $FASTQ_DIR/"$SAMPLE"_filt_paired_q20_R2_001.fastq.gz \
#               $FASTQ_DIR/"$SAMPLE"_filt_unpaired_q20_R2_001.fastq.gz \
#               AVGQUAL:20"
# ($trimmomatic) 2>&1 | tee ${OUTPUT_DIR_SAMPLE_ALIGNMENT}$sample.trimmomatic.summary.txt

#
# echo "Running bwa mem aligner: Removing sequences < mapQ20 sample- $sample"
# sudo bash /var/pipelines_ngs_test/shell/bwaAlign_exome.sh $DIR/"$sample"_filt_paired_q20_R1_001.fastq.gz   $DIR/"$sample"_filt_paired_q20_R2_001.fastq.gz   $DIR2
#
#
# echo " generating sorted bam by coordinate sample- $sample"### sorting
# java -jar /opt/picard2/picard.jar SortSam \
#           I=$DIR2/"$sample".bam  \
#           O=$DIR2/"$sample".sorted.bam  \
#           SORT_ORDER=coordinate
#
# echo "#####generating alignment stat sample- $sample"
# java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics \
#           R=/home/doc/ref/ref_genome/ucsc.hg19.fasta \
#           I=$DIR2/"$sample".sorted.bam \
#           O=$DIR2/"$sample".sorted.bam.alignmentMetrics.txt
#
# echo "#####removing duplicates sample- $sample "
# java -jar /opt/picard2/picard.jar MarkDuplicates \
#           REMOVE_DUPLICATES=true \
#           INPUT=$DIR2/"$sample".sorted.bam \
#           OUTPUT=$DIR2/"$sample".sorted.rmdups.bam \
#           METRICS_FILE=$DIR2/"$sample".sorted.rmdups.bam.metrics.txt
#
# echo "#####generating bam index sample- $sample"
# java -jar /opt/picard2/picard.jar BuildBamIndex \
#           I=$DIR2/"$sample".sorted.rmdups.bam
#
# echo "#####generating CalculateHsMetrics sample- $sample "
# # java -jar /opt/picard/picard-tools-1.134/picard.jar BedToIntervalList  I=cre_v1_design.bed O=/home/hhadmin/exome_pipeline/01_bamQC/cre_v1_design_bed.interval_list SD=/doc/ref/ref_genome/ucsc.hg19.dict
# java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
#           I=$DIR2/"$sample".sorted.rmdups.bam  \
#           O=$DIR2/"$sample".output_hs_metrics.txt \
#           R=/home/doc/ref/ref_genome/ucsc.hg19.fasta \
#           BAIT_INTERVALS= /home/hhadmin/exome_pipeline/01_bamQC/cre_v1_design_bed.interval_list \
#           TARGET_INTERVALS= /home/hhadmin/exome_pipeline/01_bamQC/cre_v1_design_bed.interval_list
#
#
# # mapping distribution among chromosome
# # /opt/samtools19/bin/samtools idxstats --threads 8 $sample.sorted.rmdups.bam > $sample.sorted.rmdups.bam.idxstats
# # /opt/python3/bin/python3 ../../03_plotIdstats.py $sample.sorted.rmdups.bam.idxstats
#
#
# # flag and map distribution
# # /opt/samtools19/bin/samtools view  --threads 8 $sample.sorted.rmdups.bam | awk '{print $2 "," $5}' > $sample.sorted.rmdups.bam.samflag.mapq
# # /opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.bam.samflag.mapq flag
# # /opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.bam.samflag.mapq mapq
#
#
# ## coverage analysis
# /opt/bedtools2/bin/bedtools bamtobed -i $DIR2/"$sample".sorted.rmdups.bam > $DIR2/"$sample".bed
# awk '$1 !~ /\_/ {print}' $DIR2/"$sample".bed > $DIR2/"$sample".filter.bed
# /opt/bedtools2/bin/bedtools coverage -a /home/hhadmin/exome_pipeline/01_bamQC/cre_v1_design.bed -b $DIR2/"$sample".filter.bed -mean > $DIR2/"$sample".CREcoverage.mean.bed
# /opt/python3/bin/python3 /home/hhadmin/exome_pipeline/01_bamQC/02_plotMeanCoverage.py $DIR2/"$sample".CREcoverage.mean.bed 500
#
#
#
#
# # /opt/python3/bin/python3 ../../05_generatePDFReport.py "$sample"
