#!/usr/bin/perl
use warnings;
use strict;

#UK10Kのvcfではindelの際 "-"を使わずに表現しているためそこをmafと同様にするためにtrimするスクリプト
open( VCF,"gunzip -c /Volumes/areca42TB/ega/file/all_sample_topdriver_region.vcf.gz|");
open(LVCF,"gunzip -c /Volumes/areca42TB/ega/file/all_sample_topdriver_region_likevcf.tsv.gz|");
open( OUT,"|gzip -c >/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle.vcf.gz");
open(LOUT,"|gzip -c >/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle_likevcf.tsv.gz");

while(<VCF>){
		if($_ =~/^#/){print OUT "$_";next;}
		chomp;
		my @line = split(/\t/,);
		if((grep{length($_) >1}($line[3],$line[4]))&&(length($line[3])==length($line[4]))){ #handle same length snp
				my ($pos,$ref,$alt)=($line[1],$line[3],$line[4]);
				while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
						$ref=substr($ref,1);
						$alt=substr($alt,1);
						$pos++;
				}
				while($ref and $alt and substr($ref,-1,1) eq substr($alt,-1,1) and $ref ne $alt){
						$ref=substr($ref,0,-1);
						$alt=substr($alt,0,-1);
				}
				print OUT "$line[0]\t$pos\t$line[2]\t$ref\t$alt\t$line[5]\t$line[6]\t$line[7]\n";
		}elsif(grep{length($_) >1}($line[3],$line[4])){ #handle same tail indel
				my ($ref,$alt)=($line[3],$line[4]);
				if(($ref=~/^$alt/)||($alt=~/^$ref/)){ #戦闘の配列が同じ(mafに変える時に自動的に直してくれるのでスルー)
						print OUT "$_\n";
				}else{
						while($ref and $alt and substr($ref,-1,1) eq substr($alt,-1,1) and $ref ne $alt and !(($ref=~/^$alt/)||($alt=~/^$ref/))){
								$ref=substr($ref,0,-1);
								$alt=substr($alt,0,-1);
						}
						print OUT "$line[0]\t$line[1]\t$line[2]\t$ref\t$alt\t$line[5]\t$line[6]\t$line[7]\n";
				}
		}else{
				print OUT "$_\n";
		}
}
close VCF;
close OUT;
		
my $header = <LVCF>;
print LOUT "$header";		
while(<LVCF>){
		chomp;
		my @line = split(/\t/,);
		if(grep{length($_) >1}($line[2],$line[3])){
				my ($pos,$ref,$alt)=@line[1..3];
				while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
						($ref,$alt) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $alt );
						$pos++;
				}
				while($ref and $alt and substr($ref,-1,1) eq substr($alt,-1,1) and $ref ne $alt){
#						if(grep{length($_)==2}($ref,$alt) and length($ref) != length($alt)){print "there are same tail indel $_\n";}
						($ref,$alt) = map{$_ = substr( $_, 0, -1 ); ( $_ ? $_ : "-" )} ( $ref, $alt );
				}
				print LOUT "$line[0]\t$pos\t$ref\t$alt\t$line[4]\t$line[5]\t$line[6]\t$line[7]\n";
				
		}else{
				print LOUT "$_\n";
		}
}
close LVCF;
close LOUT
