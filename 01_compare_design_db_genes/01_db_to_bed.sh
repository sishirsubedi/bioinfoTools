
# 1. convert database such as clinvar or cosmic in vcf file formate to appropriate format for running bedtools
# desired format is:
#                       chr      start    end     id              gene
#                       chr1    69224   69225   COSM3677745     OR4F5

DB="db_02_cosmic"
SAMPLE="CRE"
awk -v OFS=$'\t' '{ $1="chr"$1;print}' "$DB".vcf |
awk '{print $1"\t"$2"\t" $2+1 "\t" $0 }' |
grep 'GENE' |
cut -f -3,6,11 |
awk -F';' '{sub("GENE=", "", $1); print $1  }' |
awk '{split($5,a,"_"); print $1 "\t" $2 "\t" $3  "\t" $4 "\t"  a[1] }' > "$DB"_modified.bed


#run bedtools and get required information i.e. gene or postion
/opt/bedtools2/bin/bedtools intersect -a $SAMPLE.bed -b "$DB"_modified.bed -wa -wb | awk '{print $8}' > out_"$SAMPLE"_"$DB"_intersect_genes.csv
