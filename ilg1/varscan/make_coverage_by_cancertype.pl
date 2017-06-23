#!/usr/bin/perl
use strict;
use warnings;


my $bed_file="~/git/driver_genes/top_driver105.bed";
my @project=qw(hnsc ov prad thca ucec);
foreach my $project(@project){
		print "doing $project\n";
		&main($project);
}
exit;

#=================================================================
sub main ( $ ){
		my $pj=$_[0]; #project name
				my $path="/Volumes/areca42TB2/gdc/";
		my $norm_bam_dir=$path."norm_bam/$pj";
		my $response_norm = $path."varscan/$pj/response_norm.tsv";
		open(RES,$response_norm)or die "ERROR:cannot open $response_norm\n";
		
		my $colums=<RES>;chomp $colums;
		my @colum=split(/\t/,$colums);
		my($file_name,$case_id,$gender)=("","","");
		for(my $i=0;@colum>$i;$i++){
				if($colum[$i] eq "file_name"){$file_name=$i;
				}elsif($colum[$i] eq "cases_0_submitter_id"){$case_id=$i;
				}elsif($colum[$i] eq "cases_0_demographic_gender"){$gender=$i;
				}
		}
		($file_name ne "" and $case_id ne "" and $gender ne "") or die "ERROR:response norm file's colum is not normaly:$file_name,$gender,$case_id\n";
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
		mkdir $path."varscan/$pj/depth";
		my @bam_files=();
		my @pa_ids=();
		my $file_num=0;
		foreach my $pa_id(keys%info){
				if(!defined $info{$pa_id}{focal}){next;}
				push(@bam_files,"$norm_bam_dir/$info{$pa_id}{file}");
				push(@pa_ids,$pa_id);
				if(scalar(@bam_files)==50){
						$file_num++;
						my $outfile=$path."varscan/$pj/depth/out$file_num.tsv";
						open(OUT,">$outfile");
						print OUT "chr\tposition\t".join(@pa_ids)."\n";
						close OUT;
						`samtools depth -q 10 -b $bed_file @bam_files >>$outfile`;
						@bam_files=();
						@pa_ids=();
				}
		}
		$file_num++;
		my $outfile=$path."varscan/$pj/depth/out$file_num.tsv";
		open(OUT,">$outfile");
		print OUT "chr\tposition\t".join(@pa_ids)."\n";
		close OUT;
		`samtools depth -q 10 -b $bed_file @bam_files >>$outfile`;
		
		print "print coverage file of $pj\n";
		my $depth_dir=$path."varscan/$pj/depth";
		my @dpls=`ls $depth_dir|grep out`;chomp @dpls;
		my %an=(); 
		foreach my $dp_file(@dpls){
				$|=1;
				my @male=();
				open(DP,"$depth_dir/$dp_file");
				my $dpcolum=<DP>;chomp $dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				for(my $i=2;@dpcolum>$i;$i++){
						if($info{$dpcolum[$i]}{gender} eq "male"){push(@male,$i);}
				}
				while(<DP>){
						chomp;
						my @line=split(/\t/,);
						$line[0] =~s/^chr//;
						if(!defined $an{$line[0]}{$line[1]}){$an{$line[0]}{$line[1]}=0;}
						if($line[0] ne "X"){
								for(my $i=2;@line>$i;$i++){
										if($line[$i] >=8){$an{$line[0]}{$line[1]} +=2;}
								}
						}else{
								for(my $i=2;@line>$i;$i++){
										if((grep{$_ eq $i}@male)&&($line[$i] >=8)){
												$an{$line[0]}{$line[1]}++;
										}elsif((!grep{$_ eq $i}@male)&&($line[$i] >=8)){
												$an{$line[0]}{$line[1]} +=2;
										}
								}
						}
				}
				close DP;
		}
		my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);
		open(OUT,">$depth_dir/coverage_all_data_exist_patient.tsv");
		foreach my $chr(@chr){
				foreach my $posi(sort{$a<=>$b}keys %{$an{$chr}}){
						print OUT "$chr\t$posi\t$an{$chr}{$posi}\n";
				}
		}
		close OUT;
}

				
