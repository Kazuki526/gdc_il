#!/usr/bin/perl
use strict;
use warnings;

########## perl downloading.pl PROJECT ##########

if(!defined $ARGV[0]){die "ERROR::plese imput PROJECT";}
my $project = uc $ARGV[0];
my @focal_project =qw(BRCA CRC GBM HNSC KCC LGG LUAD LUSC OV PRAD THCA UCEC);
if(!grep{$_ eq $project}@focal_project){ die "ERROR::project is correct?? $project\n";}
my $pj= lc $project;
my $project_dir = $pj;
my $project_list = "\\\"TCGA-$project\\\"";
if($project eq "KCC"){$project_list ="[\\\"TCGA-KIRC\\\",\\\"TCGA-KIRP\\\",\\\"TCGA-KICH\\\"]";
}elsif($project eq "CRC"){$project_list="[\\\"TCGA-COAD\\\",\\\"TCGA-READ\\\"]";}
#sedする際にshに送られるコマンドに\が残る用に\\で\をエスケープし\"で"をエスケープしている
my ($response,$norm_manifest,$tumor_manifest,$project_json)=
	("$project_dir/$pj"."_response.tsv",
	"$project_dir/norm_bam/$pj"."norm_manifest.tsv",
	"$project_dir/tumor_bam/$pj"."tumor_manifest.tsv",
	"$project_dir/$pj"."_bam.json");


#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

#check top driver 105genes bed file exist?
my $bed="$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/mmr.json";
(-e $bed) or die "ERROR:not exist bed file:$bed\n";

#check existing bamslicing.pl
my $download_pl="$ENV{HOME}/git/gdc_il/ilg1/download/bamslicing.pl";
(-e $download_pl) or die "ERROR::download script $download_pl is not exist!!\n";
#check existing sample_bam.json
my $sample_json="$ENV{HOME}/git/gdc_il/ilg1/varscan/mmr/sample_bam.json";
(-e $sample_json) or die "ERROR::$sample_json exist??";

mkdir "$project_dir";
mkdir "$project_dir/norm_bam";
mkdir "$project_dir/tumor_bam";
`cat $sample_json|sed s/\\\"PROJECT\\\"/$project_list/ >$project_json`;
`curl --request POST --header \"Content-Type: application/json\" --data \@$project_json 'https://api.gdc.cancer.gov/files'|nkf -Lu >$response`;


#make hash for using patient focal
my $patient_list = "/Volumes/areca42TB/tcga/all_patient/patient_list.tsv";
my %patient_focal=();
open(PL,"$patient_list");
<PL>;
while(<PL>){
		chomp;
		my @line = split(/\t/,);
		$patient_focal{$line[1]}="ok";
}
close PL;

# rad response and make norm & tumor manifest
open(RES,"$response");
open(OUTN,">$norm_manifest");
open(OUTT,">$tumor_manifest");
print OUTN "id\tfilename\n";print OUTT "id\tfilename\n";
my @header = split(/\t/,<RES>);chomp @header;
my %header_num=();
for(my $i=0;@header>$i;$i++){
		if($header[$i] eq "file_name"){$header_num{'file'}=$i;
		}elsif($header[$i] eq "cases.0.submitter_id"){$header_num{'case_id'}=$i;
		}elsif($header[$i] eq "cases.0.samples.0.sample_type"){$header_num{'sample_type'}=$i;
		}elsif($header[$i] eq "id"){$header_num{'file_id'}=$i;
		}else{die "ERROR::$response have wrong header what is $header[$i]??\n";}
}
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if(! defined $patient_focal{$line[$header_num{case_id}]}){next;}
		if($line[$header_num{sample_type}] eq "Primary Tumor"){
				print OUTT "$line[$header_num{file_id}]\t$line[$header_num{file}]\n";
		}else{
				print OUTN "$line[$header_num{file_id}]\t$line[$header_num{file}]\n";
		}
}
close RES;
close OUTN;
close OUTT;


system("perl $download_pl $norm_manifest $bed $project_dir/norm_bam");
system("perl $download_pl $tumor_manifest $bed $project_dir/tumor_bam");



