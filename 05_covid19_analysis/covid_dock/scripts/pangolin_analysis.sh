#!/bin/bash
display_usage()
{
cat <<EOF >> /dev/stderr
 USAGE: $0
 OPTIONS:
 	-r Run ID
EOF
}

parse_options()
{
    IMPORT=0

		while getopts "hz:r:" opt; do

				case $opt in
					h)
					display_usage
					exit 1
					;;
	                z)
					IMPORT=$OPTARG
					;;
					r)
					RUN=$OPTARG
					;;
					:)
					echo "Option -$OPTARG requires an argument."
					;;
					\?)
					echo "Invalid option: -$OPTARG"
		   esac
    done

    if [ $IMPORT -gt 0 ] ; then
        return 0
    fi

    return 1
}

main(){

	parse_options $*
	
	### run pangolin analysis
	cd /storage/scratch/covid/container/pangolin_analysis/${RUN}
	cp /storage/scratch/covid/container/temp/Consensus.zip .
	unzip Consensus.zip
	cat Consensus/*.fasta > consensus_${RUN}_merged.fasta
	pangolin consensus_${RUN}_merged.fasta -o .
	mv lineage_report.csv lineage_report_${RUN}_v3_1.csv
	
	###combine pangolin results
	python /storage/scratch/covid/container/scripts/pangolin_analysis_combine.py

}

# ##############################################################################
# run main
# ##############################################################################
main $*


