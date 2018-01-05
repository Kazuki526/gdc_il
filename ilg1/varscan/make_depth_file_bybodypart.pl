#!/usr/bin/perl
#########this script make only tumor bam depth file #############
use strict;
use warnings;

my $bed_file="~/git/driver_genes/top_driver105.bed";
my @project=qw(brain breast lung kidney colorectal);
foreach my $project(@project){
		print "doing $project\n";
		&main($project);
}
exit;

#=================================================================
sub main ( $ ){
		my $pj=$_[0]; #project name
				my $path="/Volumes/areca42TB2/gdc/";
		my $bam_dir=$path."tumor_bam/$pj";
		my $response = $path."varscan/$pj/response_tumor.tsv";
		#read response_file and make %info{patient_id}-{file}&{gender}
		open(RES,$response)or die "ERROR:cannot open $response\n";
		my $colums=<RES>;chomp $colums;
		my @colum=split(/\t/,$colums);
		my($file_name,$case_id,$gender)=("","","");
		for(my $i=0;@colum>$i;$i++){
				if($colum[$i] eq "file_name"){$file_name=$i;
				}elsif($colum[$i] eq "cases_0_submitter_id"){$case_id=$i;
				}elsif($colum[$i] eq "cases_0_demographic_gender"){$gender=$i;
				}
		}
		($file_name ne "" and $case_id ne "" and $gender ne "") or die "ERROR:response tumor file's colum is not normaly:$file_name,$gender,$case_id\n";
		my %info=();
		while(<RES>){
				chomp;
				my @line = split(/\t/,);
				if(scalar @line != 5){
						for(my$i=0;5>=$i;$i++){
								if(!defined$line[$i] ){$line[$i] ="";}
						}
				}
				if(defined $info{$line[$case_id]}){
						$info{$line[$case_id]}{'file'}=$line[$case_id].".bam";
				}else{
						$info{$line[$case_id]}{'file'}=$line[$file_name];
						$info{$line[$case_id]}{'gender'}=$line[$gender];
						if($info{$line[$case_id]}{gender} eq ""){$info{$line[$case_id]}{gender}="female";}
				}
		}
		close RES;

		#remove errored bam by using gender_age file
		open(GA,$path."varscan/$pj/gender_age.tsv");
		<GA>; #delete colum line
		while(<GA>){
				chomp;
				my @line=split(/\t/,);
				if(!defined $info{$line[0]}){print "WARNING:there is not having bam patient $line[0]\n";}
				$info{$line[0]}{focal}=1;
		}
		close GA;
		
		print "print out $pj depth files\n";
		my @bam_files=();
		my @pa_ids=();
		my $file_num=0;
		foreach my $pa_id(keys%info){
				if(!defined $info{$pa_id}{focal}){next;}
				push(@bam_files,"$bam_dir/$info{$pa_id}{file}");
				push(@pa_ids,$pa_id);
				if(scalar(@bam_files)==50){
						$file_num++;
						my $outfile="/Volumes/areca42TB/tcga/maf_norm/$pj/depth/tout$file_num.tsv";
						open(OUT,">$outfile");
						print OUT "chr\tposition\t".join("\t",@pa_ids)."\n";
						close OUT;
						`samtools depth -q 10 -b $bed_file @bam_files >>$outfile`;
						@bam_files=();
						@pa_ids=();
				}
		}
		$file_num++;
		my $outfile="/Volumes/areca42TB/tcga/maf_norm/$pj/depth/tout$file_num.tsv";
		open(OUT,">$outfile");
		print OUT "chr\tposition\t".join("\t",@pa_ids)."\n";
		close OUT;
		`samtools depth -q 10 -b $bed_file @bam_files >>$outfile`;
}

				
