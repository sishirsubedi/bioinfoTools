#!/usr/bin/env bash

if [ $# -eq 0 ]
then
	echo "-t CELLTYPES"
  exit
fi

if test $# -gt 0
	then
	while getopts :t: opt
	do
	case $opt in
	t)
	  CELLTYPES=$OPTARG
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


for celltype in $(echo $CELLTYPES | tr "-" "\n") ; do

	echo "processing -- "$celltype"."

	data_dir="../data/${celltype}/"
	result_dir="../result/${celltype}/"
    
  	echo "step 1: run flowDensity with optimized parameters with parental gating"
  	Rscript runFlowDensity_Gated_${celltype}.R  $data_dir  $result_dir

	echo "step 2: run flowDensity with optimized parameters for comprehensive expression analysis"
  	Rscript runFlowDensity_Independent_${celltype}.R $data_dir  "${result_dir}Independent/"

	echo "step 3: PCA visualization"
  	/opt/python3/bin/python PCAanalysis_visualization.py ${celltype} "${result_dir} ${celltype}_result.csv" ${result_dir}
	
	echo "step 4: Comprehensive expression analysis visualization"
  	/opt/python3/bin/python comprehensiveExpressionAnalysis.py ${celltype} "${result_dir}Independent/ ${celltype}_result.csv" "${result_dir}/Independent/"
	
done
