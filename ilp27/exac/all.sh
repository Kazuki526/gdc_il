perl pick_driver_region.pl

perl lift_over_and_ref_check.pl


export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this 

export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$PATH
export PATH=$HOME/vep/samtools/bin:$PATH
cd ~/mskcc-vcf2maf-be943f6
perl vcf2maf_chr.pl --input-vcf /Volumes/cancer/exac/exac_nontcga_liftovered.vcf --output-maf /Volumes/cancer/exac_nontcga_topdriver.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38

perl vcf2maf_chr.pl --input-vcf /Volumes/cancer/exac/exac_nontcga_liftovered_indel.vcf --output-maf /Volumes/cancer/exac_nontcga_topdriver_indel.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38
