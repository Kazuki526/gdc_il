#!/usr/bin/perl
use warnings;
use strict;

my $pwd=`pwd`;chomp $pwd;
print "doing on $pwd\n";
if($pwd ne "/Volumes/areca42TB/exac/file"){die "ERROR:doing on wrong directory\n";}

my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my $exac_vcf="/Volumes/areca42TB/exac/file/exac_nontcga_liftovered.vcf.gz";
(-e $exac_vcf) or die "ERROR:EXAC vcf file is wrong? $exac_vcf\n";

open(VCF,"gunzip -c $exac_vcf|");
open(OUT1,"|gzip -c >exac_nontcga_liftovered_checked.vcf.gz");
open(OUT2,"|gzip -c >exac_nontcga_liftovered_checked_likevcf.tsv.gz");
print OUT2 "chr\tposi\tref\talt\tac_exac\tan_exac\n";
my ($chr,$posi,$ref19)=("","","");
my @lines=();
my @alt=();
while(<VCF>){
		if($_=~/^#/){print OUT1 $_;next;}
		chomp;
		my @line = split(/\t/,);
		if((($chr ne $line[0])||($posi != $line[1]))&&($chr ne "")){
				my $fasta = `samtools faidx $hg38 $chr:$posi-$posi`;
				my ($region,$ref)=split(/\n/,$fasta);
				if(grep{length($_) >1}($line[3],@alt)){
						if($ref19 !~ /^$ref/){print "WARNING:chr$line[0]:$posi ref_colum is not same with hg38 reference\n";}
				}else{
						if($ref19 ne $ref ){
								if(!grep{ $ref eq $_ }@alt){die "ERROR:there are not exist ref allele info at chr$chr:$posi\n";}
								my ($ac,$an)=(0,0);
								my @ac=();
								for(my $i=0;@alt>$i;$i++){
										if($lines[$i]=~/AC=(\d+);AN=(\d+)$/){
												$ac+=$1;$an=$2;push(@ac,$1);
										}else{die "ERROR: this line is error $lines[$i]\n";
										}
								}
								for(my $i=0;@alt >$i;$i++){
										my $refac=$an-$ac;
										my @ch_line=split(/\t/,$lines[$i]);
										if($alt[$i] eq $ref){
												print OUT1 "$chr\t$posi\t$ch_line[2]\t$ref\t$ref19\t$ch_line[5]\t$ch_line[6]\tAC=$refac;AN=$an\n";
												print OUT2 "$chr\t$posi\t$ref\t$ref19\t$refac\t$an\n";
										}else{
												$ch_line[3] = $ref;
												print OUT1 join("\t",@ch_line)."\n";
												print OUT2 "$chr\t$posi\t$ref\t$ch_line[4]\t$ac[$i]\t$an\n";
										}
								}
						}else{
								for(my $i=0;@alt>$i;$i++){
										my ($ac,$an);
										if($lines[$i]=~/AC=(\d+);AN=(\d+)$/){
												$ac=$1;$an=$2;
										}else{die "ERROR: this line is error $lines[$i]\n";
										}
								print OUT1 "$lines[$i]\n";
								print OUT2 "$chr\t$posi\t$ref\t$alt[$i]\t$ac\t$an\n";
								}
						}
						@lines=();
						@alt=();
				}
		}
		($chr,$posi,$ref19)=($line[0],$line[1],$line[3]);
		push(@lines,$_);
		push(@alt,$line[4]);
}
close VCF;
close OUT1;
close OUT2;

exit;
