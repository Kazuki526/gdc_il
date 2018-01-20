#!/usr/bin/perl
use strict;
use warnings;

############ usage:: perl make_coverage.pl ################

#read AF=5~% site list file AF_mid_list.tsv => make %mid_af{chr}{start}{ref}{alt}=focal is "ok"
my %mid_af=();
my $mid_af_list="/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_list.tsv";
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

open(OUTS,"|gzip -c >/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_coverage_by_patient.tsv.gz");
print OUTS "patient_id\tage\tgender\tchr\tstart\tfocal\n";
my @projects = qw(brca crc hnsc kcc lgg luad lusc ov prad thca ucec);
open(OUTC,"|gzip -c >/Volumes/areca42TB2/gdc/control_region/all_patient/coverage_all_cont.tsv.gz");
print OUTC "chr\tstart\tcancer_type\tan_white\tan_black\tan_other\n";
open(OUTX,"|gzip -c >/Volumes/areca42TB2/gdc/control/all_patient/coverage_X_male_cont.tsv.gz");
print OUTX "chr\tstart\tcancer_type\tan_white\tan_black\tan_other\n";
foreach my $project (@projects){
		print "doing $project\n";
		&make_coverage($project);
}
close OUTS;
close OUTC;
close OUTX;



exit;

#=================================================================
sub make_coverage ( $ ){
		my ($pj,$pj_dir) = ($_[0],"/Volumes/areca42TB2/gdc/control_region/$_[0]/"); #project name
		my %coverage=();
		my %coverage_xmale=();
		my @dpls=`ls $pj_dir/ndepth|grep out|grep -v tout`;
		for(my$file_num=1;@dpls >= $file_num;$file_num++){
				my %coverage_norm_focal=(); #if defined dpfocal is coverage=ok at normal sequence
				print "read $pj:depth$file_num\n";
				my($ndepth,$tdepth) = ("$pj_dir/ndepth/ndepth$file_num.tsv","$pj_dir/tdepth/tdepth$file_num.tsv");
				$|=1; #バッファのフラッシュ
		#read norm depth file
				open(DP,"$ndepth");
				my $ndpcolum=<DP>;chomp$ndpcolum;
				my @ndpcolum=split(/\t/,$ndpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						for(my $i=2;@line>$i;$i++){
								if($line[$i] >=8){$coverage_norm_focal{$line[0]}{$line[1]}{$ndpcolum[$i]}=1;}
						}
				}
				close DP;
		#read tumor depth file
				open(DP,"$tdepth");
				my $tdpcolum=<DP>;chomp$tdpcolum;
				my @tdpcolum=split(/\t/,$tdpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if($line[0] ne "X"){
								for(my $i=2;@line>$i;$i++){
										if(($line[$i] >=8)&&(defined $coverage_norm_focal{$line[0]}{$line[1]}{$tdpcolum[$i]})){
												$coverage{$line[0]}{$line[1]}{$info{$tdpcolum[$i]}{cancer_type}}{$info{$tdpcolum[$i]}{race}}+=2;
												if(defined $mid_af{"chr$line[0]"}{$line[1]}){
														print OUTS "$tdpcolum[$i]\t$info{$tdpcolum[$i]}{age}\t$info{$tdpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
												}
										}elsif(defined $mid_af{"chr$line[0]"}{$line[1]}){
												print OUTS "$tdpcolum[$i]\t$info{$tdpcolum[$i]}{age}\t$info{$tdpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
										}
								}
						}else{
								for(my $i=2;@line>$i;$i++){
										if(($info{$tdpcolum[$i]}{gender} eq "male")&&($line[$i] >=8)&&(defined $coverage_norm_focal{$line[0]}{$line[1]}{$tdpcolum[$i]})){
												$coverage_norm_focal{$tdpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}{$info{$tdpcolum[$i]}{cancer_type}}{$info{$tdpcolum[$i]}{race}}++;
												$coverage_xmale{$line[1]}{$info{$tdpcolum[$i]}{race}}++;
										}elsif(($info{$tdpcolum[$i]}{gender} ne "male")&&($line[$i] >=8)&&(defined $coverage_norm_focal{$line[0]}{$line[1]}{$tdpcolum[$i]})){
												$coverage_norm_focal{$tdpcolum[$i]}{$line[0]}{$line[1]}{tumor}="ok";
												$coverage{$line[0]}{$line[1]}{$info{$tdpcolum[$i]}{cancer_type}}{$info{$tdpcolum[$i]}{race}}+=2;
										}
										if(defined $mid_af{"chr$line[0]"}{$line[1]}){
												if(($line[$i] >=8)&&(defined $coverage_norm_focal{$line[0]}{$line[1]}{$tdpcolum[$i]})){
														print OUTS "$tdpcolum[$i]\t$info{$tdpcolum[$i]}{age}\t$info{$tdpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tok\n";
												}else{
														print OUTS "$tdpcolum[$i]\t$info{$tdpcolum[$i]}{age}\t$info{$tdpcolum[$i]}{gender}\tchr$line[0]\t$line[1]\tno\n";
												}
										}
								}
						}
				}
				close DP;
		}
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

				
