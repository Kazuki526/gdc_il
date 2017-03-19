#!/usr/bin/perl
use strict;
use warnings;

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB/tcga"){die "ERROR:doing on wrong dir\nthis perl must doing on /Volumes/areca42TB/tcga\n";}

my $bp=$ARGV[0]; #body part
my $dpdir="$pwd/maf_norm/$bp/depth";

my @files=grep{chomp;$_=~/^out\d+\.tsv$/}`ls $dpdir`;

my($topdrivers,$allexonbed)=("$ENV{HOME}/git/driver_genes/driver_genes.tsv",
							 "$ENV{HOME}/git/driver_genes/top_driverallexon.bed");
(-e $topdrivers and -e $allexonbed) or die "ERROR:driver_genes fiels are not exist";

open(IN,"$topdrivers");
my %topname=();
my $devnull=<IN>;
while(<IN>){
		chomp;
		my @line=split(/\t/,);
		if($line[1] <4){last;}
		$topname{$line[0]}="ok";
}
close IN;

open(BED,"$allexonbed");
my %topexon=();
$devnull=<BED>;
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		my ($gene,$devnull)=split(/:/,$line[3],2);
		if(!defined $topname{$gene}){next;}
		for(my $i=$line[1];$line[2] >= $i;$i++){
				$topexon{$line[0]}{$i}="$gene\t$line[1]\t$line[2]\t$line[5]"
		}
}
close BED;


open(OUT,"|gzip -c >$dpdir/exon_mean_depth_tidy.tsv.gz");
print OUT "patient_id\tchr\tgene\tstart\tend\tstrand\tmean_depth\n";
foreach my $file(@files){
		open(DP,"$dpdir/$file");
		my $colum=<DP>;chomp$colum;
		my @colum=split(/\t/,$colum);
		my ($exon,$chr)=("","");
		my @dp=();
		for(my $i=2;@colum>$i;$i++){
				$dp[$i-2]=0;
		}
		while(<DP>){
				chomp;
				my @line=split(/\t/,);
				if($line[0]=~/^chr(.+)$/){
						$chr=$1;
						if($chr eq "X"){$chr=23;}
				}else{die "ERROR:what line? $_\n";}
				if(!defined$topexon{$chr}{$line[1]}){next;}
				if($exon ne $topexon{$chr}{$line[1]}){
						if($exon ne ""){
								&print_out($chr,$exon,\@colum,\@dp);
								for(my $i=2;@colum>$i;$i++){
										$dp[$i-2]=0;
								}
						}
						$exon=$topexon{$chr}{$line[1]};
				}
				for(my $i=2;@colum>$i;$i++){
						$dp[$i-2]+=$line[$i];
				}
		}
		&print_out($chr,$exon,\@colum,\@dp);
		close DP;
}
close OUT;


#=======================================================
sub print_out{
		my($chr,$exon,$colum,$dp)=@_;
		my($gene,$start,$end,$strand)=split(/\t/,$exon);
		my $leng=$end - $start +1;
		for(my $i=2;@{$colum} > $i;$i++){
				my $mean_depth=${$dp}[$i-2] / $leng;
				print OUT "${$colum}[$i]\t$chr\t$exon\t$mean_depth\n";
		}
}
