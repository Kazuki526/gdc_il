#!/usr/bin/perl
use warnings;
use strict;


######################## perl make_depth_file.pl ###########################
my $bed_file = "/Volumes/areca42TB/GRCh38_singlefasta/control_genes_exon_with_splice_site.bed";
-e $bed_file or die "ERROR::not exist $bed_file!!\n";
my @projects = qw(brca crc gbm hnsc kcc lgg luad lusc ov prad thca ucec);
foreach my $project (@projects){
		&make_depth_file($project);
}
exit;

#===============================================================================================================================
sub make_depth_file( $ ){
		my $pj=$_[0]; #project_id
		#ASCAT data is bodypart dir
		my$tubam="$pj/tumor_bam";	#tumor bam directori
		my$nobam="$pj/norm_bam";	#norm bam directori

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
				if($_ =~ /7:(\S+)\sdownload more than 10 times so this file cannot download\?$/){
						$error_file{$1}="error";
				}
		}
		close TE;
		open(NE,"grep \"download more than 10 times so this file cannot download?\" $pj/norm_bam/result_download.txt|") or 
				die "ERROR:there is not exist normal download error list file download_errored_file.txt!!\n";
		while(<NE>){
				chomp;
				if($_ =~ /7:(\S+)\sdownload more than 10 times so this file cannot download\?$/){
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

#make depth file by each 50patients
		mkdir "$pj/tdepth";
		mkdir "$pj/ndepth";
		my @norm_bam=();
		my @tumor_bam=();
		my @patients=();
		my $file_num=0;
		foreach my $pid(keys %data){
				if((!defined $data{$pid}{tumor_bam})||(!defined $data{$pid}{norm_bam})){$data{$pid}{focal}="no";next;}else{$data{$pid}{focal}="ok";}
				push(@patients,$pid);
				if($data{$pid}{norm_bam} =~ /;/){
						push(@norm_bam,"$nobam/$pid.bam");
				}else{
						push(@norm_bam,$data{$pid}{norm_bam});
				}
				if($data{$pid}{tumor_bam} =~ /;/){
						push(@tumor_bam,"$tubam/$pid.bam");
				}else{
						push(@tumor_bam,$data{$pid}{tumor_bam});
				}
				if(scalar(@patients)==50){
						$file_num++;
						my($tdepth,$ndepth)=("$pj/tdepth/tdepth$file_num.tsv","$pj/ndepth/ndepth$file_num.tsv");
						open(OUT,">$tdepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
						open(OUT,">$ndepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
						print "make $pj:tdepth$file_num\n";
						`samtools depth -q 10 -b $bed_file @tumor_bam >>$tdepth`;
						print "make $pj:ndepth$file_num\n";
						`samtools depth -q 10 -b $bed_file @norm_bam >>$ndepth`;
						@patients=();
						@tumor_bam=();
						@norm_bam=();
				}
		}
		$file_num++;
		my($tdepth,$ndepth)=("$pj/tdepth/tdepth$file_num.tsv","$pj/ndepth/ndepth$file_num.tsv");
		open(OUT,">$tdepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
		open(OUT,">$ndepth");print OUT "chr\tposition\t",join("\t",@patients)."\n";close OUT;
		`samtools depth -q 10 -b $bed_file @tumor_bam >>$tdepth`;
		`samtools depth -q 10 -b $bed_file @norm_bam >>$ndepth`;
}
