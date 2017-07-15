#!/usr/bin/perl
use strict;
use warnings;


my %coverage=();
my %coverage_xmale=();
my %mid_af=();
my $mid_af_list="/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_list.tsv";
(-e $mid_af_list) or die "ERROR:AF_mid_list.tsv is not exist!!\n";
open(MA,$mid_af_list);
my $colum=<MA>;chomp$colum;
my @colum=split(/\t/,$colum);
my ($chrn,$posin);
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "chr"){$chrn=$i;
		}elsif($colum[$i] eq "start"){$posin=$i}
}
while(<MA>){
		chomp;
		my @line=split(/\t/,);
		$mid_af{$line[$chrn]}{$line[$posin]}="ok";
}
close MA;

open(OUTS,">/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_coverage_by_patient.tsv");
print OUTS "patient_id\tage\tgender\tchr\tstart\tfocal\n";
my @project=qw(hnsc ov prad thca ucec);
my @bodypart=qw(brain breast lung kidney colorectal);
foreach my $project(@project){
		print "doing $project\n";
		&main($project,"/Volumes/areca42TB2/gdc/varscan");
}
foreach my $bodypart (@bodypart){
		print "doing $bodypart\n";
		&main($bodypart,"/Volumes/areca42TB/tcga/maf_norm/");
}
close OUTS;

my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);
open(OUTC,">/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_all.tsv");
print OUTC "chr\tstart\tan_cancer\n";
open(OUTX,">/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_X_male.tsv");
print OUTX "chr\tstart\tan_male_cancer\n";
foreach my $chr(@chr){
		foreach my $posi(sort{$a<=>$b}keys %{$coverage{$chr}}){
				print OUTC "chr$chr\t$posi\t$coverage{$chr}{$posi}\n";
		}
		if($chr eq "X"){
				foreach my $posi (sort{$a<=>$b}keys %coverage_xmale){
						print OUTX "chr$chr\t$posi\t$coverage_xmale{$posi}\n";
				}
		}
}
close OUTC;
close OUTX;



exit;

#=================================================================
sub main ( $ ){
		my ($pj,$dp_dir) = ($_[0],"$_[1]/$_[0]/depth"); #project name
				my $path="/Volumes/areca42TB2/gdc/";

		my %patient=();
		#read gender_age and list patient
		open(GA,$path."varscan/$pj/gender_age.tsv");
		<GA>; #delete colum line
		while(<GA>){
				chomp;
				my @line=split(/\t/,$_,3);
				$patient{$line[0]}{"gender"}=($line[1] eq "" ? "male":$line[1]);
				$patient{$line[0]}{"age"}   =$line[2];
		}
		close GA;
		
		print "read depth file of $pj\n";
		#read norm depth file
		my @dpls=`ls $dp_dir|grep out|grep -v tout`;chomp @dpls;
		my %dp_focal=();
		foreach my $dp_file(@dpls){
				$|=1;
				open(DP,"$dp_dir/$dp_file");
				my $dpcolum=<DP>;chomp$dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						for(my $i=2;@line>$i;$i++){
								if($line[$i] >=8){$dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]}=1;}
						}
				}
				close DP;
		}

		#read tumor depth file
		my @dplst=`ls $dp_dir|grep tout`;chomp @dplst;
		foreach my $dp_file(@dplst){
				$|=1;
				my @male=();
				open(DP,"$dp_dir/$dp_file");
				my $dpcolum=<DP>;chomp$dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				for(my $i=2;@dpcolum>$i;$i++){
						if($patient{$dpcolum[$i]}{gender} eq "male"){push(@male,$i);}
				}
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if(!defined $coverage{$line[0]}{$line[1]}){$coverage{$line[0]}{$line[1]}=0;}
						if($line[0] ne "X"){
								for(my $i=2;@line>$i;$i++){
										if(($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
														$coverage{$line[0]}{$line[1]}+=2;
														if(defined $mid_af{"chr$line[0]"}{$line[1]}){
																print OUTS "$dpcolum[$i]\t$patient{$dpcolum[$i]}{age}\t$patient{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
														}
										}elsif(defined $mid_af{"chr$line[0]"}{$line[1]}){
														print OUTS "$dpcolum[$i]\t$patient{$dpcolum[$i]}{age}\t$patient{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
										}
								}
						}else{
								for(my $i=2;@line>$i;$i++){
										if((grep{$_ eq $i}@male)&&($line[$i] >=6)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
												$dp_focal{$dpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}++;
												$coverage_xmale{$line[1]}++;
										}elsif((!grep{$_ eq $i}@male)&&($line[$i] >=6)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
												$dp_focal{$dpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}+=2;
										}
										if(defined $mid_af{"chr$line[0]"}{$line[1]}){
												if(($line[$i] >=6)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
														print OUTS "$dpcolum[$i]\t$patient{$dpcolum[$i]}{age}\t$patient{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
												}else{
														print OUTS "$dpcolum[$i]\t$patient{$dpcolum[$i]}{age}\t$patient{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
												}
										}
								}
						}
				}
				close DP;
		}
}

				
