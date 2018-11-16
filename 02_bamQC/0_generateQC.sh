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

echo "#####generating stat " $sample
/opt/samtools19/bin/samtools flagstat --threads 5 "$sample".sort.bam > "$sample".sort.bam.flagstat

echo "#####generating mapping per chromosome " $sample
/opt/samtools19/bin/samtools idxstats --threads 5 "$sample".sort.bam > "$sample".sort.bam.idxstats
/opt/python3/bin/python3 05_plotIdstats.py "$sample".sort.bam.idxstats

/opt/samtools19/bin/samtools view  --threads 5 "$sample".sort.bam | awk '{print $2 "," $5}' > "$sample".sort.bam.samflag.mapq
/opt/python3/bin/python3 04_plotFlagMapqDistribution.py "$sample".sort.bam.samflag.mapq flag
/opt/python3/bin/python3 04_plotFlagMapqDistribution.py "$sample".sort.bam.samflag.mapq mapq

# select only 0X2 i.e. read mapped in proper pair
/opt/samtools19/bin/samtools view  --threads 5 -bf 0x2 "$sample".sort.bam | /opt/bedtools2/bin/bedtools bamtobed -i stdin > "$sample".bed

/opt/bedtools2/bin/bedtools coverage -a ../../cre_design.bed -b "$sample".bed -d > "$sample".CREcoverage.bed
/opt/python3/bin/python3 06_coveragePlot.py  "$sample".CREcoverage.bed  /home/hhadmin/exome/bamQC/coverage_genelist.csv

echo " processing coverage qc "
/opt/python3/bin/python3 07_coverageQC.py $bedcoverage_mean 300
