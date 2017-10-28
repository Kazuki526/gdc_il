#!/usr/bin/zsh
##zsh top105cancergene_maf.sh body_part patient_id
if [ $# -ne 2 ]; then
echo "needed argument is 2! $#argument was inputted" 1>&2
exit 1
fi

export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this 

cat $1/varscan_vcf/$2.indel.vcf|sed s/^chr// >$1/vcf2vep/$2.vcf
cat $1/varscan_vcf/$2.snp.vcf|sed -e '/^#/d' -e s/^chr// >>$1/vcf2vep/$2.vcf

export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$PATH
export PATH=$HOME/vep/samtools/bin:$PATH
cd ~/mskcc-vcf2maf-be943f6
perl vcf2maf_chr.pl --input-vcf /Volumes/cancer/kaz_gdc/varscan/$1/vcf2vep/$2.vcf --output-maf /Volumes/cancer/kaz_gdc/varscan/$1/maf/$2.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38

