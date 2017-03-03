#!/usr/bin/perl
use warnings;
use strict;

my $drivergene_bed_path="/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed";
if(!-e $drivergene_bed_path){die "ERROR:$drivergene_bed_path not exist!!";}
my $exac_nontcga_vcf_path="/Volumes/areca42TB/exac/file/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
if(!-e $exac_nontcga_vcf_path){die "ERROR:$exac_nontcga_vcf_path not exist!";}

#pick up only topdriver genes region
my %topdriver=();
open(BED,"$drivergene_bed_path")or die "cant open $drivergene_bed_path\n";
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0]=~ s/^chr//;
		my $symbol;
		if($line[3]=~/^([^;]+);/){$symbol=$1;}
		for(my $i=$line[1];$line[2] >=$i;$i++){
				$topdriver{$line[0]}{$i}=$symbol;
		}
}
close BED;

#write like_maf
open(VCF,"gunzip -c $exac_nontcga_vcf_path|");
open(OUT,"|gzip -c >/Volumes/areca42TB/exac/file/ExAC_nonTCGA_topdriver_like.maf.gz") or die "ERROR:output file error\n";
print OUT "chr\tstart\tend\tref\talt\tAC\tAN\tgene_symbol\tConsequence\tPolyPhen\n";
my %info=();
while(<VCF>){
		chomp;
		if($_ =~ /^##INFO=<ID=CSQ.+Format:\s(.+)">$/){
				my @info=split(/\|/,$1);
				for(my $i=0;@info > $i;$i++){
						$info{$info[$i]}=$i;
				}
				next;
		}elsif($_=~/^#/){next;}
		my @line=split(/\t/,);
		if(!defined $topdriver{$line[0]}{$line[1]}){next;} #if not drivergene region
		my ($AC,$AN,$CSQ);
		if($line[7]=~/AC=([0-9,]+);.+;AN=(\d+);.+;CSQ=(.+)$/){
				($AC,$AN,$CSQ)=($1,$2,$3);
		}else{die "what line is this\n$_\n";}
		if($line[4] !~ /,/){
				&print_out_by_variant($line[4],$AC,$AN,$CSQ,$_);
		}else{
				my @alt=split(/,/,$line[4]);
				my @AC =split(/,/,$AC);
				my @CSQ=split(/,/,$CSQ);
				for(my $t=0;@alt>$t;$t++){
						my @csq_alt;
						foreach my $vep(@CSQ){
								if($vep=~/^$alt[$t]\|/){push(@csq_alt,$vep);}
						}
						&print_out_by_variant($alt[$t],$AC[$t],$AN,join(",",@csq_alt),$_);
				}
		}
}
close VCF;
close OUT;
exit;
#====================================================
sub print_out_by_variant($ $ $ $ $){
		my($alt,$ac,$an,$csq,$line)=@_;
		my @line=split(/\t/,$line);
		my ($chr,$posi,$ref)=($line[0],$line[1],$line[3]);
		my @csq=split(/,/,$csq);
		my ($start,$end);
		my($ref_leng,$alt_leng)=(length$ref,length$alt);
		if($ref_leng != $alt_leng){
				while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
						($ref,$alt)=map{$_ =substr($_,1);($_ ? $_:"-")}($ref,$alt);
						$posi++;$ref_leng--;$alt_leng--;
				}
		}
		if($ref_leng==$alt_leng){
				($start,$end)=($posi,$posi+$alt_leng-1);
		}else{
				if($ref_leng < $alt_leng){
						($start,$end)=(($ref eq "-" ? $posi-1 : $posi),($ref eq "-" ? $posi : $posi + $ref_leng -1));
				}else{
						($start,$end)=($posi,$posi + $ref_leng - 1);
				}
		}
		for(my $i=0;@csq>$i;$i++){
				my @info=split(/\|/,$csq[$i]);
				$info[$info{Consequence}] =~ s/&/,/g;
				if($info[$info{SYMBOL}] eq $topdriver{$line[0]}{$line[1]}){
						print OUT "$chr\t$start\t$end\t$ref\t$alt\t$an\t$info[$info{SYMBOL}]\t$info[$info{Consequence}]\t$info[$info{PolyPhen}]\n";
						last;
				}elsif(scalar(@csq) -1 ==$i){
						print "WARNING: this line $alt has no meaningful variant\n$line\n";
				}
		}
}

