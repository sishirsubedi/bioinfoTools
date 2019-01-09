#!/usr/bin/env bash

# tumor="/home/hhadmin/exome_pipeline/01_bamQC/run_bcm/case_01/TCRBOA1-T-WEX.bam"
# normal="/home/hhadmin/exome_pipeline/01_bamQC/run_bcm/case_01/TCRBOA1-N-WEX.bam"
# ref="/home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa"
# sample="COV_7_COV_8"
# method="mutect"
#
#
# ####varscan somatic
tumor_pileup="/opt/samtools19/bin/samtools  mpileup  -f $ref  $tumor ";
normal_pileup="/opt/samtools19/bin/samtools  mpileup  -f $ref  $normal ";
java -jar /opt/varscan/VarScan.v2.3.9.jar somatic <($normal_pileup) <($tumor_pileup) sample_vep_output.txt

#
echo " $currentdate    INFO  -  running VEP"
/opt/vep_94/ensembl-tools-release-94/vep_94/ensembl-vep/vep \
-i /home/hhadmin/exome_pipeline/02_variantCalling/$sample/"$sample"_snp_indel_filter.vcf \
-o /home/hhadmin/exome_pipeline/02_variantCalling/$sample/"$sample"_snp_indel_filter.vep.vcf.txt \
--offline \
--dir_cache /opt/vep_94/ensembl-tools-release-94/cache \
--vcf \
--refseq \
--pick_allele \
--sift p \
--polyphen p \
--hgvs \
--symbol \
--vcf \
--pubmed \
--fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
--force_overwrite




#using samtools
# /opt/samtools19/bin/samtools view -H /home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-1/COV-1.sorted.rmdups.bam | awk -F"\t" '/@SQ/{print $2}' | cut -d":" -f2

# create chromosome bam using bamtools
/opt/bamtools/bamtools/bin/bamtools split -in  $normal -reference
/opt/bamtools/bamtools/bin/bamtools split -in  $tumor -reference

# index all chr specific bam files
# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
#   /opt/samtools19/bin/samtools index /home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-7/COV-7.sorted.rmdups.REF_chr$chr.bam
#   /opt/samtools19/bin/samtools index /home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-8/COV-8.sorted.rmdups.REF_chr$chr.bam
# done
#
for counter in 1 6 11 16 21; do

  ##check chromosome input file
    if [ -f  chr.vcf ]; then
      rm chr.vcf
    fi

  ## write to input file
    if [ $counter = "21" ]; then
      for inner_counter in 21 22 X Y;do
        chr=$inner_counter
        echo $chr >> chr.vcf
      done
    else
      for inner_counter in 0 1 2 3 4;do
        chr=$((counter + inner_counter))
        echo $chr >> chr.vcf
      done
    fi

## run gatk in parallel
cat chr.vcf | /opt/parallel/bin/parallel "java -jar /opt/GATK4/GenomeAnalysisTK.jar \
                                    -T MuTect2 \
                                    -R /home/doc/ref/ref_genome/ucsc.hg19.fasta \
                                    -I:tumor /home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-7/COV-7.sorted.rmdups.REF_chr{}.bam \
                                    -I:normal /home/hhadmin/exome_pipeline/01_bamQC/run_258_analysis/q20/COV-8/COV-8.sorted.rmdups.REF_chr{}.bam \
                                    -o chr{}_gatk.output.vcf \
                                    --min_base_quality_score 30 \
                                    --max_alt_allele_in_normal_fraction 0.00 \
                                    -L chr{}"
done


# filter gatk output files
# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
#   java -jar /opt/GATK4/GenomeAnalysisTK.jar \
#        -R /home/doc/ref/ref_genome/ucsc.hg19.fasta \
#        -T VariantsToTable \
#        -V chr"$chr"_gatk.output.vcf \
#        -F CHROM -F POS -F REF -F ALT -F FILTER \
#        -F ECNT -F HCNT -F MAX_ED -F MIN_ED -F NLOD -F RPA -F RU=CA -F STR -F TLOD \
#        -GF GT -GF AD -GF AF -GF ALT_F1R2 -GF ALT_F2R1 -GF FOXOG -GF PGT  -GF PID  -GF QSS  -GF REF_F1R2 -GF REF_F2R1 \
#        -o chr"$chr"_gatk.output.vcf.filter.txt
# done
