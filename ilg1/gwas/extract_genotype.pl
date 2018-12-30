#!/usr/bin/perl
use warnings;
use strict;

#read focal patient list
# variant_num_rank  ==>> >5%:>10, >3%:>11, >2%:>12, >1%:>14
open(PL,"/Volumes/areca42TB/tcga/all_patient/gwas_focal_patient_list.tsv");
my %focal_patient=();
my $header=<PL>;chomp $header;
my @header = split(/\t/,$header);
my %colum=();
for(my $i=0;@header>$i;$i++){$colum{$header[$i]}=$i;}
while(<PL>){
		chomp;
		my @line = split(/\t/,);
		$focal_patient{$line[$colum{"patient_id"}]}=$line[$colum{"race"}];
}
close PL;

#genotype annotation file
open(AN,"$ENV{HOME}/ascat/gw6/lib/affygw6.hg19.pfb");
<AN>;
my %annotation =();
while(<AN>){
		chomp;
		my @line = split(/\t/,);
		$annotation{$line[0]}="$line[1]:$line[2]";
}
close AN;

#cancer type ごとに
my @cancer_type=qw(brain breast colorectal hnsc kidney lung ov prad thca ucec);
foreach my $canty (@cancer_type){
#cel file name to patient_id
		print "doing $canty\n";
		open(LS,"/Volumes/areca42TB/tcga/CNA/$canty/cel/norm_cel_list.txt");
		<LS>;
		my %cel2pid=();
		while(<LS>){chomp;my@line=split(/\t/,);$cel2pid{$line[1]}=$line[2];}
		close LS;
#extract focal patient & site
		open(GE, "/Volumes/areca42TB/tcga/CNA/$canty/cel/apt_norm/birdseed.calls.txt");
		open(OUT,"|gzip -c >/Volumes/areca42TB/tcga/CNA/$canty/cel/apt_norm/extract_gentype_for_control_rare_gwas.tsv.gz");
		my@focal_colum=();
		while(<GE>){
				chomp;
				if($_=~/^#/){next;}
				my @line=split(/\t/,);
				# CEL file name line
				if($line[0] eq "probeset_id"){
						print OUT "probeset";
						for(my $i=1;@line>$i;$i++){
								if(!defined $cel2pid{$line[$i]}){die "ERROR::what CEL file?? at $canty $line[$i]\n";}
								my $pid = $cel2pid{$line[$i]};
								if(defined$focal_patient{$pid}){push(@focal_colum,$i);print OUT "\t$pid";}
						}
						print OUT "\n";
				}else{
						if(!defined $annotation{$line[0]}){next;}
						print OUT "$annotation{$line[0]}";
						foreach my $coln (@focal_colum){
								print OUT "\t$line[$coln]";
						}
						print OUT "\n";
				}
		}
		close GE;
		close OUT;
}
exit;
