#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: generateQC.sh"
	echo "-f bam file name"
	exit
fi

if test $# -gt 0
	then
	while getopts :f: opt
	do
	case $opt in
	f)
			sample=$OPTARG
			;;
	:)
		echo "Option -$OPTARG requires an argument."
		;;
	\?)
		echo "Invalid option: -$OPTARG"
	esac
	done
	shift $((OPTIND-1))
fi

echo " generating sorted bam by coordinate "### sorting
java -jar /opt/picard2/picard.jar SortSam I=$sample.bam O=$sample.sorted.bam SORT_ORDER=coordinate

echo "#####generating alignment stat $sample "
java -jar /opt/picard2/picard.jar CollectAlignmentSummaryMetrics R=/home/doc/ref/ref_genome/ucsc.hg19.fasta I=$sample.sorted.bam O=$sample.sorted.bam.alignmentMetrics.txt

echo "#####removing duplicates $sample "
java -jar /opt/picard2/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=$sample.sorted.bam OUTPUT=$sample.sorted.rmdups.bam METRICS_FILE=$sample.rmdups.metrics.txt

echo "#####generating mapq filetered bam $sample "
/opt/samtools19/bin/samtools view  --threads 10 -bf 0x2 -q 30 $sample.sorted.rmdups.bam > $sample.sorted.rmdups.filter.bam

echo "#####generating bam index $sample "
java -jar /opt/picard2/picard.jar BuildBamIndex I=$sample.sorted.rmdups.filter.bam

# mapping distribution among chromosome
/opt/samtools19/bin/samtools idxstats --threads 10 $sample.sorted.rmdups.filter.bam > $sample.sorted.rmdups.filter.bam.idxstats
/opt/python3/bin/python3 ../../03_plotIdstats.py $sample.sorted.rmdups.filter.bam.idxstats


#flag and map distribution
/opt/samtools19/bin/samtools view  --threads 10 $sample.sorted.rmdups.filter.bam | awk '{print $2 "," $5}' > $sample.sorted.rmdups.filter.bam.samflag.mapq
/opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.filter.bam.samflag.mapq flag
/opt/python3/bin/python3 ../../04_plotFlagMapqDistribution.py $sample.sorted.rmdups.filter.bam.samflag.mapq mapq


## coverage analysis
/opt/bedtools2/bin/bedtools bamtobed -i $sample.sorted.rmdups.filter.bam > "$sample".bed
awk '$1 !~ /\_/ {print}' "$sample".bed > "$sample".filter.bed
/opt/bedtools2/bin/bedtools coverage -a ../../cre_design.bed -b "$sample".filter.bed -mean > "$sample".filter.CREcoverage.mean.bed
/opt/python3/bin/python3 ../../02_plotMeanCoverage.py "$sample".filter.CREcoverage.mean.bed 500


echo "#####generating CalculateHsMetrics $sample "
# java -jar /opt/picard/picard-tools-1.134/picard.jar BedToIntervalList  I=cre_design.bed O=cre_design_bed.interval_list SD=/doc/ref/ref_genome/ucsc.hg19.dict
java -jar /opt/picard/picard-tools-1.134/picard.jar CalculateHsMetrics \
 I="$sample".sorted.bam \
 O=output_hs_metrics.txt \
 R=/home/doc/ref/ref_genome/ucsc.hg19.fasta \
 BAIT_INTERVALS= ../../cre_design_bed.interval_list \
 TARGET_INTERVALS= ../../cre_design_bed.interval_list

/opt/python3/bin/python3 ../../05_generatePDFReport.py "$sample"
