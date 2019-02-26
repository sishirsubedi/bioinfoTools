#!/usr/bin/env bash

tumor="/home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-1/COV-1.sorted.rmdups.bam"
normal="/home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-2/COV-2.sorted.rmdups.bam"
ref="/home/doc/ref/ref_genome/ucsc.hg19.fasta"
outdir="/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/ebcall/"
# rundir="/home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/ebcall/EBCall-master/"
sh ebCall_v2.sh $tumor $normal $outdir /home/hhadmin/exome_pipeline/02_variantCalling/COV_1_COV_2/ebcall/gen_ref.txt
