#!/usr/bin/perl
use warnings;
use strict;

######################## perl varscan.pl PROJECT ###########################

my $pj=$ARGV[0];#project_id
#ASCAT data is bodypart dir
my %bodypart=("hnsc" =>"hnsc","ov"=>"ov","prad"=>"prad","thca"=>"thca","ucec"=>"ucec",
		"brca"=>"breast","crc"=>"colorectal","kcc"=>"kidney","luad"=>"lung","lusc"=>"lung","gbm"=>"brain","lgg"=>"brain");
my$tubam="$pj/tumor_bam";	#tumor bam directori
my$nobam="$pj/norm_bam";	#norm bam directori

#check reference of bam are exist?
my $ref="/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";

#check pj ascat file exist ?
my $ascat_file="/Volumes/areca42TB/tcga/CNA/$bodypart{$pj}/cel/annotate_ascat.tsv.gz";
(-e $ascat_file)or die "ERROR:not exist ascat file:$ascat_file\n";

#make hash for using patient focal
my $patient_list = "/Volumes/areca42TB/tcga/all_patient/patient_list.tsv";
my %patient_focal=();
open(PL,"$patient_list") or die "ERROR:cannot open $patient_list\n";
<PL>;
while(<PL>){
		chomp;
		my @line = split(/\t/,);
		$patient_focal{$line[1]}="ok";
}
close PL;
#check download errored file
my %error_file=();
open(TE,"grep \"download more than 10 times so this file cannot download?\" $pj/tumor_bam/result_download.txt|") or 
die "ERROR:there is not exist tumor download error list file download_errored_file.txt!!\n";
while(<TE>){
		chomp;
		if($_ =~ /\d+:(\S+)\sdownload more than 10 times so this file cannot download\?$/){
				$error_file{$1}="error";
		}
}
close TE;
open(NE,"grep \"download more than 10 times so this file cannot download?\" $pj/norm_bam/result_download.txt|") or 
die "ERROR:there is not exist normal download error list file download_errored_file.txt!!\n";
while(<NE>){
		chomp;
		if($_ =~ /\d+:(\S+)\sdownload more than 10 times so this file cannot download\?$/){
				$error_file{$1}="error";
		}
}
close NE;

my %data=();
my $response = "$pj/$pj"."_response.tsv";
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
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if((!defined $patient_focal{$line[$header_num{case_id}]})||(defined$error_file{$line[$header_num{file}]})){next;}
		if($line[$header_num{sample_type}] eq "Primary Tumor"){
				if(defined $data{$line[$header_num{case_id}]}{tumor_bam}){
						$data{$line[$header_num{case_id}]}{tumor_bam}.=";$tubam/$line[$header_num{file}]";
				}else{
						$data{$line[$header_num{case_id}]}{tumor_bam}="$tubam/$line[$header_num{file}]";
				}
		}
		else{
				if(defined $data{$line[$header_num{case_id}]}{norm_bam}){
						$data{$line[$header_num{case_id}]}{norm_bam}.=";$nobam/$line[$header_num{file}]";
				}else{
						$data{$line[$header_num{case_id}]}{norm_bam}="$nobam/$line[$header_num{file}]";
				}
		}
}
close RES;

#get purity from ascat data
open(ASCAT,"gunzip -c $ascat_file|")or die "ERROR:cannot open $ascat_file\n";
<ASCAT>;#0:gene_symbol 1:patient_id 2:nmajor 3:nminor 4:chr 5:cna_start 6:cna_end 7:start_rate  8:end_rate 9:gene_start 10:gene_end 11:ploidy 12:purity
my $patient_id="";
while(<ASCAT>){
		chomp;
		my @line=split(/\t/,);
		if(($line[1] eq $patient_id)||(!defined$patient_focal{$line[1]})){next;}
		$patient_id=$line[1];
		$data{$patient_id}{'purity'}=$line[12];
}
close ASCAT;

#merge bam files to 1bam file when 1patient have many bam
print "merging bam files\n";
foreach my $pid(keys %data){
		if((!defined $data{$pid}{tumor_bam})||(!defined $data{$pid}{norm_bam})){$data{$pid}{focal}="no";next;}else{$data{$pid}{focal}="ok";}
		if($data{$pid}{norm_bam} =~ /;/){
				my @bam=split(/;/,$data{$pid}{norm_bam});
				`samtools merge -f $nobam/$pid.bam @bam`;
				`samtools index $nobam/$pid.bam`;
				$data{$pid}{'file_norm'}="$nobam/$pid.bam";
		}else{
				`samtools index $data{$pid}{norm_bam}`;
				$data{$pid}{'file_norm'}=$data{$pid}{norm_bam};
		}
		if($data{$pid}{tumor_bam} =~ /;/){
				my @bam=split(/;/,$data{$pid}{tumor_bam});
				`samtools merge -f $tubam/$pid.bam @bam`;
				`samtools index $tubam/$pid.bam`;
				$data{$pid}{'file_tumor'}="$tubam/$pid.bam";
		}else{
				`samtools index $data{$pid}{tumor_bam}`;
				$data{$pid}{'file_tumor'}=$data{$pid}{tumor_bam};
		}
}

mkdir "$pj/vcf";
my$num=0;
foreach my $pid(keys %data){
		$num++;
		if((!defined$data{$pid}{focal})||($data{$pid}{focal} eq "no")){print "$pid cannot download both tumor & normal files\n";next;}
		print "$num:varscan $pid\n";
		my $mpile = "samtools mpileup -q 10 -f $ref $data{$pid}{file_norm} $data{$pid}{file_tumor}";
		`zsh -c \"varscan somatic <\($mpile\) $pj/vcf/$pid --tumor-purity $data{$pid}{purity} --p-value 0.1 --output-vcf 1 --mpileup 1\"`;
}
exit;

