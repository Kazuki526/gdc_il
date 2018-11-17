#!/usr/bin/perl
use strict;
use warnings;

### useage: perl gene_list2bed_json.pl splice:NUM


my $mmr_list="$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/mmr.tsv";
#splice site 前後数塩基をbedに含ませたい時用
my $splice_length = 0;
if(defined $ARGV[0]){
		if($ARGV[0] =~ /^splice:(\d+)$/){$splice_length = $1;
		}else{die "ERROR::what input in ARGV = $ARGV[0] ?? this means splice region length\n";}
}
-e $mmr_list or die "ERROR::enter the MMR genes list? or file name are correct  $mmr_list?\n";
my %mmr_genes=();
open(IN,"$mmr_list") or die "ERROR::cannot open $mmr_list\n";
<IN>;#header
while(<IN>){
		chomp;
		my@line=split(/\t/,);
		$mmr_genes{$line[0]}=$line[1];
}
close IN;

my $gff="/Volumes/areca42TB/Homo_sapiens.GRCh38.90.gff3.gz";
-e $gff or die "ERROR::not exist GRCh38 gff3 file\n";

open(OUT,">$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/mmr_exon_with_splice_site.bed");
my ($focal,$gene)=(0,"");
open(GFF,"gunzip -c $gff|") or die "ERROR::cannot open $gff\n";
while(<GFF>){
		if($_ =~ /^#/){$focal=0;$gene="";next;}
		if($focal==0){
				chomp;
				my @line=split(/\t/,);
				if($line[2] eq "biological_region"){next;}
				if(($line[2] ne "gene")||($line[8] !~/;biotype=protein_coding;/)){$focal=2;next;}
				if($line[8] =~/;Name=([^;]+);/){$gene=$1;}else{die "ERROR::GFF have no gene name?\n$_\n";}
				if(defined $mmr_genes{$gene}){$focal=1;
				}else{$focal=2;}
		}elsif($focal==2){next;
		}else{
				my @line=split(/\t/,);
				if(($line[2] eq "exon")&&($line[8] =~/^Parent=transcript:(ENST[^;]+);/)){
						if($1 ne $mmr_genes{$gene}){next;}
						$line[3]--; #bedではここの長さがgffと違うため
						$line[3] -= $splice_length;$line[4] += $splice_length; #splice regon の数塩基分をbedに含ませるため
						print OUT "chr$line[0]\t$line[3]\t$line[4]\t$gene:$mmr_genes{$gene}\t0\t$line[6]\n";
				}
		}
}
close GFF;

close OUT;
exit;
