#!/usr/bin/perl
use warnings;
use strict;

my $single_tsv="/Volumes/areca42TB/GRCh38_singlefasta/GRCh38_single_proteincoding.tsv";
-e $single_tsv or die "ERROR::not existGRCh38_single_proteincoding.tsv\n";

my $gff="/Volumes/areca42TB/Homo_sapiens.GRCh38.90.gff3.gz";
-e $gff or die "ERROR::not exist GRCh38 gff3 file\n";

my $driver_genes_file="$ENV{HOME}/git/innanlab/driver_genes.tsv";
-e $driver_genes_file or die "ERROR::$driver_genes_file is not exist\n";
open(DG,"$driver_genes_file");
my%driver_genes=();
<DG>;
while(<DG>){
		chomp;
		my @line=split(/\t/,);
		$driver_genes{$line[0]}="CG";
}
close DG;


#pick single copy gene_symbol and ENST
open (IN,"$single_tsv") or die "ERROR::cannot open $single_tsv\n";
my %singene=();
my %sinenst=();
while(<IN>){
		chomp;
		my @line=split(/\t/,);
		$line[0]=~s/^chr//;
		$line[1]++;#######################################################################################
		$singene{$line[3]}{chr}=$line[0];
		$sinenst{$line[4]}.="$line[1]-$line[2];";
		if(defined $singene{$line[3]}{enst}){
				if($singene{$line[3]}{enst} !~ /$line[4]/){
						$singene{$line[3]}{enst}.="$line[4];";
				}
		}else{$singene{$line[3]}{enst}="$line[4];";}
}
close IN;

my %gene=();
my %enst=();
my %strand=();
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
				$gene{$line[0]}{$gene}{start}=$line[3];
				$gene{$line[0]}{$gene}{end}=$line[4];
				$gene{$line[0]}{$gene}{strand}=$line[6];
				$focal=1;
		}elsif($focal==2){next;
		}else{
				my @line=split(/\t/,);
				if($line[2] eq "mRNA"){
						if($line[8] =~/^ID=transcript:([^;]+);/){
								$gene{$line[0]}{$gene}{enst}.="$1;";
						}else{die "ERROR::GFF this line of mRNA has no ID??\n$_\n";
						}
				}elsif($line[2] eq "CDS"){
						if($line[8] =~ /;Parent=transcript:([^;]+);/){
								$enst{$1}.="$line[3]-$line[4];";
								$strand{$1}=$line[6];
						}else{die "ERROR::GFF line of CDS has no transcript ID:\n$_\n;";
						}
				}
		}
}
close GFF;

foreach my $gene (keys %singene){
		my $chr=$singene{$gene}{chr};
		if(!defined $gene{$chr}{$gene}{enst}){die "ERROR::GFF has no data about $gene\n";
		}else{
				my @enst=sort(split(/;/,$gene{$chr}{$gene}{enst}));
				my @sinenst=sort(split(/;/,$singene{$gene}{enst}));
				if(scalar(@enst) != scalar(@sinenst)){
						$gene{$chr}{$gene}{focal}="no";
				}else{
						my$gene_focal=0;
						for(my$i=0;@enst>$i;$i++){
								if($enst[$i] ne $sinenst[$i]){$gene_focal++;}
						}
						if(($gene_focal==0)&&(&enst_check($singene{$gene}{enst}))){
								$gene{$chr}{$gene}{focal}="ok";
						}else{
								$gene{$chr}{$gene}{focal}="no";
						}
				}
		}
}

foreach my $chr (sort keys %gene){
		my ($gene_before,$end)=("",0);
		foreach my $gene (sort {$gene{$chr}{$a}{start} <=> $gene{$chr}{$b}{start}}keys %{$gene{$chr}}){
				if($gene{$chr}{$gene}{start} <= $end ){$gene{$chr}{$gene}{focal}="no";$gene{$chr}{$gene_before}{focal}="no";
				}elsif(!defined$gene{$chr}{$gene}{focal}){$gene{$chr}{$gene}{focal}="no";
				}
				if($gene{$chr}{$gene}{end} > $end){$end = $gene{$chr}{$gene}{end};$gene_before=$gene;}
		}
}

open(OUT,">/Volumes/areca42TB/GRCh38_singlefasta/single_nonoverlap.tsv");
foreach my $chr (sort keys %gene){
		if($chr eq "Y"){next;}
		foreach my $gene (sort {$gene{$chr}{$a}{start} <=> $gene{$chr}{$b}{start}}keys %{$gene{$chr}}){
				if((!defined $driver_genes{$gene})&&($gene{$chr}{$gene}{focal} eq "ok")){
						print OUT "$gene\t$chr\t$gene{$chr}{$gene}{start}\t$gene{$chr}{$gene}{end}\n";
				}
		}
}

#######################################################################################################################
sub enst_check(){
		my @ensts = split(/;/,$_[0]);
		my$focal=0;
		foreach my $enst(@ensts){
				$sinenst{$enst}=&stop_codon_add($strand{$enst},$sinenst{$enst});
				if($enst{$enst} ne $sinenst{$enst}){$focal++;}
		}
		if($focal > 0){return 0;
		}else{return 1;}
}



########################################################################################################################
sub stop_codon_add(){
		my ($strand,$st_end)=@_;
		my @start=();my @end=();
		foreach(split(/;/,$st_end)){
				my @st_end=split(/-/,);
				push(@start,$st_end[0]);
				push(@end,$st_end[1]);
		}
		@start=sort{$a <=> $b}@start;
		@end = sort{$a <=> $b}@end;
		if($strand eq "+"){$end[$#end] += 3;
		}elsif($strand eq "-"){$start[0]-=3;
		}else{die "ERROR::given strand is wrong? $strand\n";}
		my $out="";
		for(my $i=0;@start>$i;$i++){
				$out.="$start[$i]-$end[$i];";
		}
		return($out);
}





