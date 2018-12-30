#!/usr/bin/perl
use warnings;
use strict;

#cancer type ごとに
open(BA,"gunzip -c /Volumes/areca42TB/tcga/CNA/brain/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(BR,"gunzip -c /Volumes/areca42TB/tcga/CNA/breast/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(CR,"gunzip -c /Volumes/areca42TB/tcga/CNA/colorectal/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(HN,"gunzip -c /Volumes/areca42TB/tcga/CNA/hnsc/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(KI,"gunzip -c /Volumes/areca42TB/tcga/CNA/kidney/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(LU,"gunzip -c /Volumes/areca42TB/tcga/CNA/lung/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(OV,"gunzip -c /Volumes/areca42TB/tcga/CNA/ov/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(PR,"gunzip -c /Volumes/areca42TB/tcga/CNA/prad/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(TH,"gunzip -c /Volumes/areca42TB/tcga/CNA/thca/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(UC,"gunzip -c /Volumes/areca42TB/tcga/CNA/ucec/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz|");
open(OUT,"|gzip -c >/Volumes/areca42TB/tcga/array_genotype/gentype_for_control_rare_gwas.tsv.gz");
while(<BA>){
		chomp;
		print OUT $_;
		my $line = <BR>; &print_line($line);
		$line = <CR>; &print_line($line);
		$line = <HN>; &print_line($line);
		$line = <KI>; &print_line($line);
		$line = <LU>; &print_line($line);
		$line = <OV>; &print_line($line);
		$line = <PR>; &print_line($line);
		$line = <TH>; &print_line($line);
		$line = <UC>; &print_line($line);
		print OUT "\n";
}
close BA;
close BR;
close CR;
close HN;
close KI;
close LU;
close OV;
close PR;
close TH;
close UC;
close OUT;


sub print_line( $ ){
		my $line = $_[0];
		chomp $line;
		if($line =~/^[\dXYMT]+:\d+\t/){print OUT "\t$'";}
		elsif($line =~ /^probeset\t/){print OUT "\t$'";}
}
