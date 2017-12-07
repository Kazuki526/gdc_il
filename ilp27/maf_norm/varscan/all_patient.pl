#!/usr/local/perl
use strict;
use warnings;

#perl all_patient.pl body_part(or project)

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc/varscan"){die "ERROR: doing on wrong dir!!\n";}

my $bp = $ARGV[0]; #body part

mkdir "$bp/vcf2vep";
mkdir "$bp/maf";

(-e "$bp/gender_age.tsv") or die "ERROR: there is not gender_age.tsv!!\n";
open(GA,"$bp/gender_age.tsv") or die "ERROR: cant open gender_age.tsv!!\n";
<GA>;
while(<GA>){
		chomp;
		my @line=split(/\t/,);
		(-e "$bp/varscan_vcf/$line[0].indel.vcf" and -e "$bp/varscan_vcf/$line[0].snp.vcf") or die "WARNING:there are not rsynced varscan vcf!!\n";
		system("sh ~/git/gdc_il/ilp27/maf_norm/varscan/varscan_vcf2maf.sh $bp $line[0]");
}
close GA;
exit;
		
