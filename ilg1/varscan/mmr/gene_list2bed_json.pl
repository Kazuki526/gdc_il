#!/usr/bin/perl
use strict;
use warnings;

### useage: perl gene_list2bed_json.pl

my $mmr_list="$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/mmr.tsv";
-e $mmr_list or die "ERROR::enter the MMR genes list? or file name are correct  $mmr_list?\n";
my %cont_genes=();
open(IN,"$mmr_list") or die "ERROR::cannot open $mmr_list\n";
<IN>;#header
while(<IN>){
		chomp;
		my@line=split(/\t/,);
		$cont_genes{$line[0]}=$line[1];
}
close IN;

my $gff="/Volumes/areca42TB/Homo_sapiens.GRCh38.90.gff3.gz";
-e $gff or die "ERROR::not exist GRCh38 gff3 file\n";

open(OUT,">$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/mmr.json");
print OUT "{\n\t\"regions\":[\n";
my ($focal,$gene)=(0,"");
open(GFF,"gunzip -c $gff|") or die "ERROR::cannot open $gff\n";
my $times=0;
while(<GFF>){
		if($_ =~ /^#/){$focal=0;$gene="";next;}
		if($focal==0){
				chomp;
				my @line=split(/\t/,);
				if($line[2] eq "biological_region"){next;}
				if(($line[2] ne "gene")||($line[8] !~/;biotype=protein_coding;/)){$focal=2;next;}
				if($line[8] =~/;Name=([^;]+);/){$gene=$1;}else{die "ERROR::GFF have no gene name?\n$_\n";}
				if(defined $cont_genes{$gene}){$focal=1;
				}else{$focal=2;}
		}elsif($focal==2){next;
		}else{
				my @line=split(/\t/,);
				if(($line[2] eq "exon")&&($line[8] =~/^Parent=transcript:(ENST[^;]+);/)){
						if($1 ne $cont_genes{$gene}){next;}
						$line[3]--;
						if($times==0){print OUT "\t\t\"chr$line[0]:$line[3]-$line[4]\"";
						}else{print OUT ",\n\t\t\"chr$line[0]:$line[3]-$line[4]\"";
						}
						$times++;
				}
		}
}
close GFF;

print OUT "\n\t]\n}";
close OUT;
exit;
