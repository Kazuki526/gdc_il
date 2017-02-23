#!/usr/bin/perl
use warnings;
use strict;

#doing on : perl push_genotype_maf.pl BODYPART
my $bp=$ARGV[0];

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(10);
my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $vcfdir="$bp/samtools";
my $mafdir="$bp/maf";

open(LIST,"$bp/list_of_perfect_maf.tsv") or die "ERROR list_of_perfect_maf.tsv exist?\n";
my$devnull=<LIST>;
my %patients=();
while(<LIST>){
		chomp;
		my @line=split(/\t/,);
		$patients{$line[0]}=$line[1];
}
close LIST;

mkdir "./$bp/genotyped_maf";

foreach my $pa (keys %patients){
		if($patients{$pa} eq "no"){next;}
#fork and returns the pid for child
		my $pid=$pm->start and next;
		print "doing $pa\n";
		open(VCF,"$vcfdir/$pa.vep.vcf") or die "$pa.vcf is not exist\n";
		my %vcf=();
		while(<VCF>){
				if($_=~/^#/){next;}
				my @line=split(/\t/,);
				my($depth,$ref_count,$alt_count);
				if($line[7] =~ /DP=(\d+);.+;DP4=([^;]+);/){
						$depth=$1;
						my @count=split(/,/,$2);
						if(scalar(@count) !=4){die "ERROR $pa.vcf:$_\nthis line have error DP4 so edit the script\n";}
						$ref_count = $count[0] + $count[1];
						$alt_count = $count[2] + $count[3];
				}
				my($allele1,$allele2,$ref,$alt,$start,$end);
				if($line[9] =~ /^(\d)\/(\d):/){
						my($genotype1,$genotype2)=($1,$2);
						if(($genotype1==2)||($genotype2==2)){die "ERROR $pa.vcf:$_\nthis line have error genotype so edit the script\n";}
						my $posi = $line[1];
						$ref=$line[3];
						if($line[4] =~ /^([ATGC]+),/){$alt=$1;
						}else{$alt = $line[4];}
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
						if($genotype1==0){$allele1=$ref;
						}else{$allele1=$alt;}
						$allele2=$alt;

				}
				$vcf{$line[0]}{$start}{$end}{dp}=$depth;
				$vcf{$line[0]}{$start}{$end}{rc}=$ref_count;
				$vcf{$line[0]}{$start}{$end}{ac}=$alt_count;
				$vcf{$line[0]}{$start}{$end}{a1}=$allele1;
				$vcf{$line[0]}{$start}{$end}{a2}=$allele2;
				$vcf{$line[0]}{$start}{$end}{ref}=$ref;
		}
		close VCF;
		open(MAF,"$mafdir/$pa.maf") or die "$pa.maf is not exist\n";
		open(OUT,">$bp/genotyped_maf/$pa.maf");
		$devnull=<MAF>;
		my $colum_name=<MAF>;
		my @colum=split(/\t/,$colum_name);
		if(($colum[11] ne "Tumor_Seq_Allele1")||($colum[12] ne "Tumor_Seq_Allele2")||($colum[39] ne "t_depth")||($colum[40] ne "t_ref_count")||($colum[41] ne "t_alt_count")){
				die "ERROR $pa.maf colum name is error\n";
		}
		print OUT "$devnull$colum_name";
		while(<MAF>){
				chomp;
				my @line=split(/\t/,);
				if($line[10] ne $vcf{$line[4]}{$line[5]}{$line[6]}{ref}){print "ERROR push_genotype_maf.pl line54-59 is mistake in $pa!!\n$_\n";next;}
				$line[11] = $vcf{$line[4]}{$line[5]}{$line[6]}{a1};
				$line[12] = $vcf{$line[4]}{$line[5]}{$line[6]}{a2};
				$line[39] = $vcf{$line[4]}{$line[5]}{$line[6]}{dp};
				$line[40] = $vcf{$line[4]}{$line[5]}{$line[6]}{rc};
				$line[41] = $vcf{$line[4]}{$line[5]}{$line[6]}{ac};
				my $out;
				if(scalar(@colum) > scalar(@line)){
						$out=join("\t",@line);
						for(my $i=0;scalar(@colum)-scalar(@line)>$i;$i++){
								$out.="\t";
						}
						$out.="\n";
				}else{
						$out=join("\t",@line)."\n";
				}
				print OUT "$out";
		}
		close MAF;
		close OUT;

		$pm->finish;#terminates the chld.process
}
exit;
				
