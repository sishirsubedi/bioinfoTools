#!/usr/bin/env bash
export SHELL=/usr/bin/bash
#
# if [ $# -eq 0 ]
# then
# 	echo "Usage: generateQC.sh"
# 	echo "-f bam file name"
# 	exit
# fi
#
# if test $# -gt 0
# 	then
# 	while getopts :f: opt
# 	do
# 	case $opt in
# 	f)
# 			sample=$OPTARG
# 			;;
# 	:)
# 		echo "Option -$OPTARG requires an argument."
# 		;;
# 	\?)
# 		echo "Invalid option: -$OPTARG"
# 	esac
# 	done
# 	shift $((OPTIND-1))
# fi
sample="COV-7"
DIR="/home/nextseq/181029_NS500761_0258_AH2JH2BGX9/out1"
DIR2="/home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20"
# echo "Running Trimmomatic: Removing sequences < Q20 sample- $sample"
# sudo java -jar /opt/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -threads 15 \
#               $DIR/"$sample"_S5_R1_001.fastq.gz \
#               $DIR/"$sample"_S5_R2_001.fastq.gz \
#               $DIR/"$sample"_filt_paired_q20_S5_R1_001.fastq.gz \
#               $DIR/"$sample"_filt_unpaired_q20_S5_R1_001.fastq.gz \
#               $DIR/"$sample"_filt_paired_q20_S5_R2_001.fastq.gz \
#               $DIR/"$sample"_filt_unpaired_q20_S5_R2_001.fastq.gz \
#               AVGQUAL:20

echo "Running bwa mem aligner: Removing sequences < mapQ20 sample- $sample"
sudo bash /var/pipelines_ngs_test/shell/bwaAlign_exome.sh $DIR/"$sample"_filt_paired_q20_S5_R1_001.fastq.gz   $DIR/"$sample"_filt_paired_q20_S5_R2_001.fastq.gz   $DIR2/"$sample"


echo " generating sorted bam by coordinate sample- $sample"### sorting
java -jar /opt/picard2/picard.jar SortSam \
          I=$DIR2/"$sample"/"$sample".bam  \
          O=$DIR2/"$sample"/"$sample".sorted.bam  \
          SORT_ORDER=coordinate

echo "#####generating alignment stat sample- $sample"
java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics \
          R=/home/doc/ref/ref_genome/ucsc.hg19.fasta \
          I=$DIR2/"$sample"/"$sample".sorted.bam \
          O=$DIR2/"$sample"/"$sample".sorted.bam.alignmentMetrics.txt

echo "#####removing duplicates sample- $sample "
java -jar /opt/picard2/picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=true \
          INPUT=$DIR2/"$sample"/"$sample".sorted.bam \
          OUTPUT=$DIR2/"$sample"/"$sample".sorted.rmdups.bam \
          METRICS_FILE=$DIR2/"$sample"/"$sample".sorted.rmdups.bam.metrics.txt

echo "#####generating bam index sample- $sample"
java -jar /opt/picard2/picard.jar BuildBamIndex \
          I=$DIR2/"$sample"/"$sample".sorted.rmdups.bam

echo "#####generating CalculateHsMetrics sample- $sample "
# java -jar /opt/picard/picard-tools-1.134/picard.jar BedToIntervalList  I=cre_design.bed O=cre_design_bed.interval_list SD=/doc/ref/ref_genome/ucsc.hg19.dict
java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
          I=$DIR2/"$sample"/"$sample".sorted.rmdups.bam  \
          O=$DIR2/"$sample"/"$sample".output_hs_metrics.txt \
          R=/home/doc/ref/ref_genome/ucsc.hg19.fasta \
          BAIT_INTERVALS= /home/hhadmin/exome_pipeline/01_bamQC/cre_design_bed.interval_list \
          TARGET_INTERVALS= /home/hhadmin/exome_pipeline/01_bamQC/cre_design_bed.interval_list


# mapping distribution among chromosome
# /opt/samtools19/bin/samtools idxstats --threads 8 $sample.sorted.rmdups.bam > $sample.sorted.rmdups.bam.idxstats
# /opt/python3/bin/python3 ../../03_plotIdstats.py $sample.sorted.rmdups.bam.idxstats


# flag and map distribution
# /opt/samtools19/bin/samtools view  --threads 8 $sample.sorted.rmdups.bam | awk '{print $2 "," $5}' > $sample.sorted.rmdups.bam.samflag.mapq
# /opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.bam.samflag.mapq flag
# /opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.bam.samflag.mapq mapq


## coverage analysis
/opt/bedtools2/bin/bedtools bamtobed -i $DIR2/"$sample"/"$sample".sorted.rmdups.bam > $DIR2/"$sample"/"$sample".bed
awk '$1 !~ /\_/ {print}' $DIR2/"$sample"/"$sample".bed > $DIR2/"$sample"/"$sample".filter.bed
/opt/bedtools2/bin/bedtools coverage -a /home/hhadmin/exome_pipeline/01_bamQC/cre_design.bed -b $DIR2/"$sample"/"$sample".filter.bed -mean > $DIR2/"$sample"/"$sample".CREcoverage.mean.bed
/opt/python3/bin/python3 /home/hhadmin/exome_pipeline/01_bamQC/02_plotMeanCoverage.py $DIR2/"$sample"/"$sample".CREcoverage.mean.bed 500




# /opt/python3/bin/python3 ../../05_generatePDFReport.py "$sample"
