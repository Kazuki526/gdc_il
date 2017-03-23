#!/usr/bin/perl
use warnings;
use strict;

my @chr=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X");
my %vcf=();
my $header="";
my @sample=();
my @vcf_main_colums=();
foreach my $chr (@chr){
		print "doing chr$chr\n";
		open(VCF,"/Volumes/cancer/1kg/topdrivers_liftovered/chr$chr.vcf");
		while(<VCF>){
				chomp;
				if($_=~/^##/){
						if($chr eq "1"){$header.="$_\n";
						}elsif($_=~/^##bcftools/){$header.="$_\n";}
				}elsif($_=~/^#CHROM/){
						my @line=split(/\t/,);
						@sample=splice(@line,9,$#line);
						if($chr eq "X"){$header.=join("\t",@line[0..8]);}
				}else{
						my @line=split(/\t/,);
						if($line[4]=~/[^ATGC,]/){next;}
						push(@vcf_main_colums,join("\t",@line[0..8]));
						for(my $i=0;@sample>$i;$i++){
								my $colum=$i+9;
								if($line[$colum] ne "0|0"){
										$vcf{$sample[$i]}{"$line[0]:$line[1]"}=$line[$colum];
								}
						}
				}
		}
		close VCF;
}

foreach my $sample (@sample){
		print "printing $sample\n";
		open(OUT1,">/Volumes/cancer/1kg/by_sample/$sample"."_allele1.vcf");
		open(OUT2,">/Volumes/cancer/1kg/by_sample/$sample"."_allele2.vcf");
		print OUT1 "$header\t$sample"."_allele1\n";
		print OUT2 "$header\t$sample"."_allele2\n";
		foreach my $vcf_main(@vcf_main_colums){
				my @vcf_main=split(/\t/,$vcf_main);
				my $key="$vcf_main[0]:$vcf_main[1]";
				if(defined $vcf{$sample}{$key}){
						my @gt=split(/\|/,$vcf{$sample}{$key});
						if($vcf_main[4] =~ /,/){
								my @alt =split(/,/,$vcf_main[4]);
								if($gt[0] !=0){
										my $gti=$gt[0] -1;
										$vcf_main[4]=$alt[$gti];
										print OUT1 join("\t",@vcf_main)."\t$gt[0]\n";
								}
								if(scalar(@gt) ==1){next;
								}elsif($gt[1] !=0){
										my $gti=$gt[1] -1;
										$vcf_main[4]=$alt[$gti];
										print OUT2 join("\t",@vcf_main)."\t$gt[1]\n";
								}
						}else{
								if($gt[0] ==1){
										print OUT1 "$vcf_main\t$gt[0]\n";
								}
								if(scalar(@gt) ==1){next;
								}elsif($gt[1] ==1){
										print OUT2 "$vcf_main\t$gt[1]\n";
								}
						}
				}
		}
		close OUT1;
		close OUT2;
		my ($alle1,$alle2)=("$sample"."_allele1","$sample"."_allele2");
		system("zsh /Volumes/cancer/1kg/vcf2maf.sh $alle1");
		system("zsh /Volumes/cancer/1kg/vcf2maf.sh $alle2");
}

