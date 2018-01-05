#!/usr/bin/perl
use strict;
use warnings;

############ usage:: perl make_coverage.pl ################

my %mid_af=();
#read AF=5~% site list file AF_mid_list.tsv => make %mid_af{chr}{start}{ref}{alt}=focal is "ok"
my $mid_af_list="/Volumes/areca42TB/tcga/all_patient/AF_mid_list.tsv";
(-e $mid_af_list) or die "ERROR:AF_mid_list.tsv is not exist!!\n";
open(MA,$mid_af_list);
my $colum=<MA>;chomp$colum;
my @colum=split(/\t/,$colum);
my ($chrn,$posin,$refn,$altn); #colum position of "chr" & "start" & "ref" & "alt"
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "chr"){$chrn=$i;
		}elsif($colum[$i] eq "start"){$posin=$i;
		}elsif($colum[$i] eq "ref"  ){$refn =$i;
		}elsif($colum[$i] eq "alt"  ){$altn =$i;
		}
}
while(<MA>){
		chomp;
		my @line=split(/\t/,);
		$mid_af{$line[$chrn]}{$line[$posin]}{$line[$refn]}{$line[$altn]}="ok";
}
close MA;

#all patient info => %info
my %info=();
open(INFO,"/Volumes/areca42TB/tcga/all_patient/all_patient_racegenderage.tsv");
$colum=<INFO>;chomp $colum;
@colum=split(/\t/,$colum);
my($patidn,$cantyn,$racen,$gendern,$agen);
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "patient_id"){$patidn=$i;
		}elsif($colum[$i] eq "cancer_type"){$cantyn=$i;
		}elsif($colum[$i] eq "race"    ){$racen =$i;
		}elsif($colum[$i] eq "gender"  ){$gendern =$i;
		}elsif($colum[$i] eq "age"     ){$agen=$i;
		}
}
while(<INFO>){
		chomp;
		my @line = split(/\t/,);
		$info{$line[$patidn]}{cancer_type} = ($line[$cantyn ] eq "" ? "NA" : $line[$cantyn ]);
		$info{$line[$patidn]}{race       } = ($line[$racen  ] eq "" ? "NA" : $line[$racen  ]);
		$info{$line[$patidn]}{gender     } = ($line[$gendern] eq "" ? "NA" : $line[$gendern]);
		$info{$line[$patidn]}{age        } = ($line[$agen   ] eq "" ? "NA" : $line[$agen   ]);
}

open(OUTS,">/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_coverage_by_patient.tsv");
print OUTS "patient_id\tage\tgender\tchr\tstart\tfocal\n";
my @project=qw(hnsc ov prad thca ucec);
my @bodypart=qw(brain breast lung kidney colorectal);
open(OUTC,">/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_all.tsv");
print OUTC "chr\tstart\tcancer_type\tan_white\tan_black\tan_other\n";
open(OUTX,">/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_X_male.tsv");
print OUTX "chr\tstart\tcancer_type\tan_white\tan_black\tan_other\n";
foreach my $project(@project){
		print "doing $project\n";
		&main($project,"/Volumes/areca42TB2/gdc/varscan");
}
foreach my $bodypart (@bodypart){
		print "doing $bodypart\n";
		&main($bodypart,"/Volumes/areca42TB/tcga/maf_norm/");
}
close OUTS;
close OUTC;
close OUTX;



exit;

#=================================================================
sub main ( $ ){
		my ($pj,$dp_dir) = ($_[0],"$_[1]/$_[0]/depth"); #project name
				my $path="/Volumes/areca42TB2/gdc/";
		print "read depth file of $pj\n";
		my %coverage=();
		my %coverage_xmale=();
		#read norm depth file
		my @dpls=`ls $dp_dir|grep out|grep -v tout`;chomp @dpls;
		my %dp_focal=();
		foreach my $dp_file(@dpls){
				$|=1; #バッファのフラッシュ
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
				$|=1; #バッファのフラッシュ
				my @male=();
				open(DP,"$dp_dir/$dp_file");
				my $dpcolum=<DP>;chomp$dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				for(my $i=2;@dpcolum>$i;$i++){
						if($info{$dpcolum[$i]}{gender} eq "male"){push(@male,$i);}
				}
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if(!defined $coverage{$line[0]}{$line[1]}){$coverage{$line[0]}{$line[1]}=0;}
						if($line[0] ne "X"){
								for(my $i=2;@line>$i;$i++){
										if(($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
														$coverage{$line[0]}{$line[1]}{$info{$dpcolum[$i]}{cancer_type}}{$info{$dpcolum[$i]}{race}}+=2;
														if(defined $mid_af{"chr$line[0]"}{$line[1]}){
																print OUTS "$dpcolum[$i]\t$info{$dpcolum[$i]}{age}\t$info{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
														}
										}elsif(defined $mid_af{"chr$line[0]"}{$line[1]}){
														print OUTS "$dpcolum[$i]\t$info{$dpcolum[$i]}{age}\t$info{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
										}
								}
						}else{
								for(my $i=2;@line>$i;$i++){
										if((grep{$_ eq $i}@male)&&($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
												$dp_focal{$dpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}{$info{$dpcolum[$i]}{cancer_type}}{$info{$dpcolum[$i]}{race}}++;
												$coverage_xmale{$line[1]}{$info{$dpcolum[$i]}{race}}++;
										}elsif((!grep{$_ eq $i}@male)&&($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
												$dp_focal{$dpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}{$info{$dpcolum[$i]}{cancer_type}}{$info{$dpcolum[$i]}{race}}+=2;
										}
										if(defined $mid_af{"chr$line[0]"}{$line[1]}){
												if(($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
														print OUTS "$dpcolum[$i]\t$info{$dpcolum[$i]}{age}\t$info{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
												}else{
														print OUTS "$dpcolum[$i]\t$info{$dpcolum[$i]}{age}\t$info{$dpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
												}
										}
								}
						}
				}
				close DP;
				my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);
				foreach my $chr(@chr){
						foreach my $posi(sort{$a<=>$b}keys %{$coverage{$chr}}){
								foreach my $cancer_type (keys%{$coverage{$chr}{$posi}}){
										if(!defined($coverage{$chr}{$posi}{$cancer_type}{white})){$coverage{$chr}{$posi}{$cancer_type}{white}=0;}
										if(!defined($coverage{$chr}{$posi}{$cancer_type}{black})){$coverage{$chr}{$posi}{$cancer_type}{black}=0;}
										if(!defined($coverage{$chr}{$posi}{$cancer_type}{other})){$coverage{$chr}{$posi}{$cancer_type}{other}=0;}
										print OUTC "chr$chr\t$posi\t$cancer_type\t$coverage{$chr}{$posi}{$cancer_type}{white}\t$coverage{$chr}{$posi}{$cancer_type}{black}\t$coverage{$chr}{$posi}{$cancer_type}{other}\n";
								}
						}
						if($chr eq "X"){
								foreach my $posi(sort{$a<=>$b}keys %{$coverage{$chr}}){
										foreach my $cancer_type (keys%{$coverage{$chr}{$posi}}){
												if(!defined($coverage{$chr}{$posi}{$cancer_type}{white})){$coverage{$chr}{$posi}{$cancer_type}{white}=0;}
												if(!defined($coverage{$chr}{$posi}{$cancer_type}{black})){$coverage{$chr}{$posi}{$cancer_type}{black}=0;}
												if(!defined($coverage{$chr}{$posi}{$cancer_type}{other})){$coverage{$chr}{$posi}{$cancer_type}{other}=0;}
												print OUTX "chr$chr\t$posi\t$cancer_type\t$coverage{$chr}{$posi}{$cancer_type}{white}\t$coverage{$chr}{$posi}{$cancer_type}{black}\t$coverage{$chr}{$posi}{$cancer_type}{other}\n";
										}
								}
						}
				}
		}
}

				
