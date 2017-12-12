#!/usr/bin/perl
use strict;
use warnings;

########## perl downloading.pl PROJECT ##########

my $project = uc $ARGV[0];
my @focal_project =qw(BRCA CRC GBM HNSC KCC LGG LUAD LUSC OV PRAD THCA UCEC);
if(!grep{$_ eq $project}@focal_project){ die "ERROR::project is correct?? $project\n";}
my $pj= lc $project;
my $project_list = "\\\"TCGA-$project\\\"";
if($project eq "KCC"){$project_list ="[\\\"TCGA-KIRC\\\",\\\"TCGA-KIRP\\\",\\\"TCGA-KICH\\\"]";
}elsif($project eq "CRC"){$project_list="[\\\"TCGA-COAD\\\",\\\"TCGA-READ\\\"]";}
#sedする際にshに送られるコマンドに\が残る用に\\で\をエスケープし\"で"をエスケープしている
my ($response,$project_json)=
	("$pj/$pj"."_response.tsv","$pj/$pj"."_bam.json");


#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

#check existing sample_bam.json
my $sample_json="$ENV{HOME}/git/gdc_il/ilg1/varscan/control_gene_region/sample_bam.json";
(-e $sample_json) or die "ERROR::$sample_json exist??";

mkdir "$pj";
mkdir "$pj/norm_bam";
mkdir "$pj/tumor_bam";
`cat $sample_json|sed s/\\\"PROJECT\\\"/$project_list/ >$project_json`;
`curl --request POST --header \"Content-Type: application/json\" --data \@$project_json 'https://gdc-api.nci.nih.gov/files'|nkf -Lu >$response`;


#make hash for using patient focal
my $patient_list = "/Volumes/areca42TB/tcga/all_patient/patient_list.tsv";
my %patient_focal=();
open(PL,"$patient_list") or die "ERROR::cannot open $patient_list\n";
<PL>;
while(<PL>){
		chomp;
		my @line = split(/\t/,);
		$patient_focal{$line[1]}="ok";
}
close PL;

# rad response and make norm & tumor manifest
open(RES,"$response");
my @header = split(/\t/,<RES>);chomp @header;
my %header_num=();
for(my $i=0;@header>$i;$i++){
		if($header[$i] eq "file_name"){$header_num{'file'}=$i;
		}elsif($header[$i] eq "cases.0.submitter_id"){$header_num{'case_id'}=$i;
		}elsif($header[$i] eq "cases.0.samples.0.sample_type"){$header_num{'sample_type'}=$i;
		}elsif($header[$i] eq "id"){$header_num{'file_id'}=$i;
		}else{die "ERROR::$response have wrong header what is $header[$i]??\n";}
}
my @group=qw(1 2 3 4 5 6 7);
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if(! defined $patient_focal{$line[$header_num{case_id}]}){next;}
		if($line[$header_num{sample_type}] eq "Primary Tumor"){
				my @files=map{if(!-e "$pj$_"){die "ERROR:dir $pj$_ is not exist!!\n"}else{"$pj$_/tumor_bam/$line[$header_num{file_id}]"}@group;
				`samtools merge $pj/tumor_bam/$line[$header_num{file_id}] @files`;
		}else{
				my @files=map{if(!-e "$pj$_"){die "ERROR:dir $pj$_ is not exist!!\n"}else{"$pj$_/norm_bam/$line[$header_num{file_id}]"}@group;
				`samtools merge $pj/norm_bam/$line[$header_num{file_id}] @files`;
		}
}
close RES;

open(OUTN,">$pj/norm_bam/result_download.txt");
open(OUTT,">$pj/tumor_bam/result_download.txt");
my%error_file=();
foreach my $group (@group){
		open(TE,"grep \"download more than 10 times so this file cannot download?\" $pj$group/tumor_bam/result_download.txt|") or 
				die "ERROR:there is not exist tumor download error list file download_errored_file.txt!!\n";
		while(<TE>){
				chomp;
				if($_ =~ /\d+:(\S+)\sdownload more than 10 times so this file cannot download?$/){
						$error_file{tumor_bam}{$1}++;
				}
		}
		close TE;
		open(NE,"grep \"download more than 10 times so this file cannot download?\" $pj$group/norm_bam/result_download.txt|") or 
				die "ERROR:there is not exist normal download error list file download_errored_file.txt!!\n";
		while(<NE>){
				chomp;
				if($_ =~ /\d+:(\S+)\sdownload more than 10 times so this file cannot download?$/){
						$error_file{norm_bam}{$1}=++;
				}
		}
		close NE;
}

foreach my$nt(keys %error_file){
		foreach my$file(keys%{$erro_file{$nt}}){
				if($nt="norm_bam"){
						print OUTN "$error_file{$nt}{$file}:$file download more than 10 times so this file cannot download?\n";
				}else{
						print OUTT "$error_file{$nt}{$file}:$file download more than 10 times so this file cannot download?\n";
				}
		}
}
close OUTT;
close OUTN;



