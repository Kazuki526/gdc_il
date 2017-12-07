#!/usr/bin/perl
use warnings;
use strict;

my $single_tsv="/Volumes/areca42TB/GRCh38_singlefasta/GRCh38_single_proteincoding.tsv";
-e $single_tsv or die "ERROR::not existGRCh38_single_proteincoding.tsv\n";

my $gff="/Volumes/areca42TB/Homo_sapiens.GRCh38.84.gff3.gz";
-e $gff or die "ERROR::not exist GRCh38 gff3 file\n";


#pick single copy gene_symbol and ENST
open (IN,"$single_tsv") or die "ERROR::cannot open $single_tsv\n";
my %singene=();
my %sinenst=();
while(<IN>){
		chomp;
		my @line=split(/\t/,);
		$line[0]=~s/^chr//;
		$singene{$line[3]}{chr}=$line[0];
		$singene{$line[3]}{enst}.="$line[4];";
		$sinenst{$line[4]}.="$line[1]-$line[2];";
}
close IN;

my %gene=();
my %enst=();
my ($focal,$gene)=(0,"");
open(GFF,"gunzip -c $gff|") or die "ERROR::cannot open $gff\n";
while(<GFF>){
		if($_ =~ /^#/){$focal=0;$gene="";next;}
		if($focal==0){
				chomp;
				my @line=split(/\t/,);
				if(($line[2] ne "gene")||($line[8] !~/;biotype=protein_coding;/)){$focal=2;next;}
				if($line[8] =~/;Name=([^;]+);/){$gene=$1;}else{die "ERROR::GFF have no gene name?\n$_\n";}
				$gene{$gene}{chr}=$line[0];
				$gene{$gene}{start}=$line[3];
				$gene{$gene}{end}=$line[4];
		}elsif($focal==1){next;
		}else{
				if($line[2] eq "mRNA"){
						if($line[8] =~/^ID=transcript:([^;]);/){
								$gene{$gene}{enst}.="$1;";
						}else{die "ERROR::GFF this line of mRNA has no ID??\n$_\n";
						}
				}elsif($line[2] eq "CDS"){
						if($line[8] =~ /;Parent=transcript:([^;]+);/){
								$enst{$1}.="$line[3]-$line[4];";
						}else{die "ERROR::GFF line of CDS has no transcript ID:\n$_\n;";
						}
				}
		}
}
close GFF;

