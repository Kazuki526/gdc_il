#doing on focal dir
export OUT=`pwd`

perl ~/git/gdc_il/ilp27/exac/control_region/pick_region.pl single_nonoverlap.tsv

perl ~/git/gdc_il/ilp27/exac/control_region/lift_over_and_ref_check.pl


export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this 

export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$PATH
export PATH=$HOME/vep/samtools/bin:$PATH
cd ~/mskcc-vcf2maf-be943f6
perl vcf2maf_chr.pl --input-vcf $OUT/exac_nontcga_liftovered.vcf --output-maf $OUT/exac_nontcga_control_region.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38

perl vcf2maf_chr.pl --input-vcf $OUT/exac_nontcga_liftovered_indel.vcf --output-maf $OUT/exac_nontcga_control_region_indel.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38
