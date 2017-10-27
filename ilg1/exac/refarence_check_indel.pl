#!/usr/bin/perl
use warnings;
use strict;

my $pwd=`pwd`;chomp $pwd;
print "doing on $pwd\n";
if($pwd ne "/Volumes/areca42TB/exac/file"){die "ERROR:doing on wrong directory\n";}

my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my $exac_vcf="/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_indel.vcf.gz";
(-e $exac_vcf) or die "ERROR:EXAC vcf file is wrong? $exac_vcf\n";

open(VCF,"gunzip -c $exac_vcf|");
open(OUT1,"|gzip -c >exac_nontcga_liftovered_checked_indel.vcf.gz");
open(OUT2,"|gzip -c >exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz");
print OUT2 "chr\tposi\tref\talt\tac_exac\tan_exac\n";
my ($chr,$posi,$ref19)=("","","");
my @lines=();
my @alt=();
while(<VCF>){
		if($_=~/^#/){print OUT1 $_;next;}
		chomp;
		my @line = split(/\t/,);
		my ($ac,$an);
		if($line[7] =~/^AC=(\d+);AN=(\d+)$/){
				$ac=$1;
				$an=$2;
		}else{die "ERROR: what line!? $_\n";
		}
		print OUT1 "$_\n";
		print OUT2 "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$ac\t$an\n";
}
close VCF;
close OUT1;
close OUT2;


my $exac_vcfh="/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_indel_nothandled.vcf.gz";
(-e $exac_vcfh) or die "ERROR:EXAC vcf file is wrong? $exac_vcf\n";

open(VCFH,"gunzip -c $exac_vcfh|");
open(OUT,"|gzip -c >exac_nontcga_liftovered_checked_indel_nothandled.vcf.gz");
my ($chr,$posi,$ref19)=("","","");
my @lines=();
my @alt=();
while(<VCFH>){
		if($_=~/^#/){print OUT1 $_;next;}
		chomp;
		my @line = split(/\t/,);
		my ($ac,$an);
		if($line[7] =~/^AC=(\d+);AN=(\d+)$/){
				$ac=$1;
				$an=$2;
		}else{die "ERROR: what line!? $_\n";
		}
		print OUT1 "$_\n";
}
close VCFH;
close OUT;

exit;
