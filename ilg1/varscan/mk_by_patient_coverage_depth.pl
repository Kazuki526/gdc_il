#!/usr/bin/perl
use strict;
use warnings;

############ usage:: perl make_coverage.pl ################
#read bed file to make focal site is oncogene or TSG
my %site_role = ();   #$site_role{CHR}{POSI} = ROLE
my %driver_gene = (); #$driver_gene{GENE} = ROLE
open(DG,"$ENV{HOME}/git/driver_genes/driver_genes.tsv") or die "ERROR::cannot open driver_genes.tsv\n";
<DG>; #delete header line
while(<DG>){
		chomp;
		my @line = split(/\t/,);
		if($line[1] < 4){last;}
		if($line[2] eq "oncogene/TSG"){$line[2] = "TSG";}
		$driver_gene{$line[0]} = $line[2];
}
close DG;
my $all_exon_bp = 0;
open(BED,"$ENV{HOME}/git/driver_genes/top_driverallexon.bed") or die "ERROR::connot open top_driverallexon.bed\n";
<BED>; #delete header line
while(<BED>){
		chomp;
		my @line = split(/\t/,);
		my $gene_symbol ="";
		if($line[3] =~ /^([^:]+):ENST\d+:\d+$/){$gene_symbol = $1;}
	#	if(!defined $driver_gene{$gene_symbol}){next;}
		for(my $posi = $line[1] +1; $line[2] >= $posi; $posi++){
				$site_role{$line[0]}{$posi} = $driver_gene{$gene_symbol};
				$all_exon_bp++;
		}
}
close BED;


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

open(OUT,"|gzip -c >/Volumes/areca42TB2/gdc/varscan/all_patient/by_patient_coverage_depth.tsv.gz");
print OUT "patient_id\tcancer_type\trace\tage\tstage\trole\tnorm_mean_depth\ttumor_mean_depth\tcalled_bp\n";
my @project=qw(ov hnsc prad thca ucec);
my @bodypart=qw(brain breast lung kidney colorectal);
foreach my $project(@project){
		print "doing $project\n";
		&main($project,"/Volumes/areca42TB2/gdc/varscan");
}
foreach my $bodypart (@bodypart){
		print "doing $bodypart\n";
		&main($bodypart,"/Volumes/areca42TB/tcga/maf_norm/");
}
close OUT;



exit;

#=================================================================
sub main ( $ ){
		my ($pj,$dp_dir) = ($_[0],"$_[1]/$_[0]/depth"); #project name
				my $path="/Volumes/areca42TB2/gdc/";
		print "read depth file of $pj\n";
		my %mean_depth = ();
		my %coverage_bp = (); #coverage_by{patient_id}{cancer_gene role} = coverage bp
		#read norm depth file
		my @dpls_normal=`ls $dp_dir|grep out|grep -v tout`;chomp @dpls_normal;
		my %dp_focal=(); #if defined dpfocal is coverage=ok at normal sequence
		foreach my $dp_file(@dpls_normal){
				print "read $pj:$dp_file\n";
				$|=1; #バッファのフラッシュ
				open(DP,"$dp_dir/$dp_file") or die "ERROR::cannot open $dp_dir/$dp_file\n";
				my $dpcolum=<DP>;chomp$dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if(!defined $site_role{$line[0]}{$line[1]}){next;}
						for(my $i=2;@line>$i;$i++){
								if($line[$i] >=8){$dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]}=1;}
								$mean_depth{$dpcolum[$i]}{norm}{$site_role{$line[0]}{$line[1]}} += $line[$i];
						}
				}
				close DP;
		}
		#read tumor depth file
		my @dpls_tumor=`ls $dp_dir|grep tout`;chomp @dpls_tumor;
		foreach my $dp_file(@dpls_tumor){
				print "read $pj:$dp_file\n";
				$|=1; #バッファのフラッシュ
				open(DP,"$dp_dir/$dp_file")or die "ERROR::cannot open $dp_dir/$dp_file\n";
				my $dpcolum=<DP>;chomp$dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if(!defined $site_role{$line[0]}{$line[1]}){next;}
						for(my $i=2;@line>$i;$i++){
								if(($line[$i] >=8)&&(defined $dp_focal{$line[0]}{$line[1]}{$dpcolum[$i]})){
										$coverage_bp{$dpcolum[$i]}{$site_role{$line[0]}{$line[1]}}++;
								}
								$mean_depth{$dpcolum[$i]}{tumor}{$site_role{$line[0]}{$line[1]}} += $line[$i];
						}
				}
				close DP;
		}
		foreach my $patient( keys %coverage_bp){
				$mean_depth{$patient}{norm}{oncogene}  /= $all_exon_bp;
				$mean_depth{$patient}{tumor}{oncogene} /= $all_exon_bp;
				$mean_depth{$patient}{norm}{TSG}       /= $all_exon_bp;
				$mean_depth{$patient}{tumor}{TSG}      /= $all_exon_bp;
				if(!defined $all_info{$patient}){$all_info{$patient} = "not_reported";}
				print OUT "$patient\t$info{$patient}{cancer_type}\t$info{$patient}{race}\t$info{$patient}{age}\t";
				print OUT "$all_info{$patient}\toncogene\t$mean_depth{$patient}{norm}{oncogene}\t";
				print OUT "$mean_depth{$patient}{tumor}{oncogene}\t$coverage_bp{$patient}{oncogene}\n";
				print OUT "$patient\t$info{$patient}{cancer_type}\t$info{$patient}{race}\t$info{$patient}{age}\t";
				print OUT "$all_info{$patient}\tTSG\t$mean_depth{$patient}{norm}{TSG}\t";
				print OUT "$mean_depth{$patient}{tumor}{TSG}\t$coverage_bp{$patient}{TSG}\n";
		}
}

				
