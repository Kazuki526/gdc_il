#!/usr/bin/perl
use strict;
use warnings;

############ usage:: perl make_coverage.pl ################
#all patient(bam called patients) info => %info
my %info=();
open(INFO,"/Volumes/areca42TB/tcga/all_patient/all_patient_racegenderage.tsv");
my $colum=<INFO>;chomp $colum;
my @colum=split(/\t/,$colum);
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
		my $cancer_type;
		if($line[$cantyn] eq "KIRC"| $line[$cantyn] eq "KIRP" | $line[$cantyn] eq "KICH"){
				$cancer_type ="KCC";
		}elsif($line[$cantyn] eq "COAD" | $line[$cantyn] eq "READ"){
				$cancer_type = "CRC";
		}else{$cancer_type = $line[$cantyn];}
		$info{$line[$patidn]}{cancer_type} = ($cancer_type    eq "" ? "NA" : $cancer_type   );
		$info{$line[$patidn]}{race       } = ($line[$racen  ] eq "" ? "NA" : $line[$racen  ]);
		$info{$line[$patidn]}{gender     } = ($line[$gendern] eq "" ? "NA" : $line[$gendern]);
		$info{$line[$patidn]}{age        } = ($line[$agen   ] eq "" ? "NA" : $line[$agen   ]);
}
close INFO;

#all patient info (all patient) => %all_info
my %all_info = ();
open(AINFO,"$ENV{HOME}/git/all_patient/all_patient_response.tsv") or die "ERROR::cannot open all_patient_info.tsv\n";
$colum =<AINFO>;chomp $colum;
@colum = split(/\t/,$colum);
my ($paidn,$stagen);
for(my $i=0;@colum > $i;$i++){
		if($colum[$i] eq "submitter_id"){$paidn = $i;
		}elsif($colum[$i] eq "diagnoses.0.tumor_stage"){$stagen = $i;}
}
while(<AINFO>){
		chomp;
		my @line = split(/\t/,);
		$all_info{$line[$paidn]} = $line[$stagen];
}
close AINFO;

mkdir "all_patient";
open(OUT,"|gzip -c >all_patient/by_patient_coverage_depth_mmr.tsv.gz") or die "ERROR::print out file have error\n";
print OUT "patient_id\tcancer_type\trace\tage\tstage\tnorm_mean_depth\ttumor_mean_depth\tcalled_bp\n";
my @projects = qw(brca crc gbm hnsc kcc lgg luad lusc ov prad thca ucec);
foreach my $project(@projects){
		print "doing $project\n";
		&main($project);
}
close OUT;



exit;

#=================================================================
sub main ( $ ){
		my ($pj,$pj_dir) = ($_[0],"$_[0]"); #project name
		print "read depth file of $pj\n";
		#read norm depth file
		my @dpls=`ls $pj_dir/ndepth|grep ndepth`;
		for(my$file_num=1;@dpls >= $file_num;$file_num++){
				my %mean_depth = ();
				my %coverage_bp = ();
				my %dp_focal=(); #if defined dpfocal is coverage=ok at normal sequence
				print "read $pj:depth$file_num\n";
				my($ndepth,$tdepth) = ("$pj_dir/ndepth/ndepth$file_num.tsv","$pj_dir/tdepth/tdepth$file_num.tsv");
				$|=1; #バッファのフラッシュ
		#read norm depth file
				open(DP,"$ndepth") or die "ERROR::cannot open $ndepth\n";
				my $ndpcolum=<DP>;chomp$ndpcolum;
				my @ndpcolum=split(/\t/,$ndpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						for(my $i=2;@line>$i;$i++){
								if($line[$i] >=8){$dp_focal{$line[0]}{$line[1]}{$ndpcolum[$i]}=1;}
								$mean_depth{$ndpcolum[$i]}{norm} += $line[$i];
						}
				}
				close DP;
		#read tumor depth file
				open(DP,"$tdepth") or die "ERROR::cannot open $tdepth\n";
				my $tdpcolum=<DP>;chomp$tdpcolum;
				my @tdpcolum=split(/\t/,$tdpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						for(my $i=2;@line>$i;$i++){
								if(($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$tdpcolum[$i]})){
										$coverage_bp{$tdpcolum[$i]}++;
								}
								$mean_depth{$tdpcolum[$i]}{tumor} += $line[$i];
						}
				}
				close DP;
				foreach my $patient( keys %coverage_bp){
						print OUT "$patient\t$info{$patient}{cancer_type}\t$info{$patient}{race}\t$info{$patient}{age}\t";
						print OUT "$all_info{$patient}\t$mean_depth{$patient}{norm}\t$mean_depth{$patient}{tumor}\t$coverage_bp{$patient}\n";
				}
		}
}

				
