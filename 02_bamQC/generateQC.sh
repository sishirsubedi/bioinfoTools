#!/usr/bin/env bash
export SHELL=/usr/bin/bash

if [ $# -eq 0 ]
then
	echo "Usage: generateQC.sh"
	echo "-f bam file name"
	echo "-b bedcoverage"
  echo "-m bedcoverage mean"
  echo "-c chromosome map"
  echo "-f flagstats mapq"
	exit
fi

if test $# -gt 0
	then
	while getopts :f:b:m:c:f: opt
	do
	case $opt in
		f)
			sample=$OPTARG
			;;
  b)
		bedcoverage=$OPTARG
		;;
    m)
      bedcoverage_mean=$OPTARG
      ;;
    c)
  		chrmap=$OPTARG
  		;;
      f)
    		flagstat_mapq=$OPTARG
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

###generate fastqc metrics
/opt/fastqc/FastQC/fastqc $sample.bam

echo " processing coverage qc "
/opt/python3/bin/python3 07_coverageQC.py $bedcoverage_mean 300
echo " processing bed coverage plot for chromosome  "
/opt/python3/bin/python3 06_coveragePlot.py $bedcoverage chr1
echo " processing coverage histogram  "
/opt/python3/bin/python3 05_compareChrmap.py $chrmap
echo " processing map quality histogram  "
/opt/python3/bin/python3 04_compareMapq.py $flagstat_mapq mapq
