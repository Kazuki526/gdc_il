#!/usr/bin/zsh
##zsh top105cancergene_maf.sh body_part bamfile_name patient_id
if [ $# -ne 3 ]; then
echo "needed argument is 3! $#argument was inputted" 1>&2
exit 1
fi

export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this 

bam=$1/$1_norm_cancergenes_bam/$2
samtools mpileup -uf /Volumes/cancer/GRCh38.d1.vd1.fa -l ~/git/gdc_il/ilp27/maf_norm/top_driver105.bed $bam|bcftools call -O v -v -c|vcf-annotate -f +/d=5/|grep -e 'PASS' -e '^#'|sed s/chr// >$1/samtools/$3.vcf


export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$PATH
export PATH=$HOME/vep/samtools/bin:$PATH
cd ~/mskcc-vcf2maf-be943f6
perl vcf2maf.pl --input-vcf /Volumes/cancer/kaz_gdc/$1/samtools/$3.vcf --output-maf /Volumes/cancer/kaz_gdc/$1/maf/$3.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --ncbi-build GRCh38

