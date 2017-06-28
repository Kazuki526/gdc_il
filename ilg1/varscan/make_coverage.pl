#!/usr/bin/perl
use strict;
use warnings;


my %coverage=();
my %coverage_xmale=();
my %mid_af_coverage=();
my $mid_af_list="AF_mid_list.tsv";
(-e $mid_af_list) or die "ERROR:AF_mid_list.tsv is not exist!!\n";
open(MA,$mid_af_list);
my @colum=split(/\t/,<MA>);
my ($chrn,$posin);
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "chr"){$chn=$i;
		}elsif($colum[$i] eq "start"){$posin=$i}
}
while(<MA>){
		chomp;
		my @line=split(/\t/,);
		$mid_af_coverage{$line[$chrn]}{$line[$posin]}{'focal'}="ok";
}
close MA;

my @project=qw(hnsc ov prad thca ucec);
my @bodypart=qw(brain breast lung kidney colorectal);
foreach my $project(@project,"/Volumes/areca42TB/gdc/varscan"){
		print "doing $project\n";
		&main($project);
}
exit;

#=================================================================
sub main ( $ ){
		my ($pj,$dp_dir) = ($_[0],"$_[1]/$_[0]/depth"); #project name
				my $path="/Volumes/areca42TB2/gdc/";

		my %patient=();
		#remove errored bam by using gender_age file
		open(GA,$path."varscan/$pj/gender_age.tsv");
		<GA>; #delete colum line
		while(<GA>){
				chomp;
				my @line=split(/\t/,);
				$patient{$line[0]}{"gender"}=$line[1];
				$patient{$line[0]}{"age"}   =$line[2];
		}
		close GA;
		
		print "print coverage file of $pj\n";
		my $dp_dir=$path."varscan/$pj/depth";
		my @dpls=`ls $dp_dir|grep out|grep -v tout`;chomp @dpls;
		foreach my $dp_file(@dpls){
				$|=1;
				my @male=();
				open(DP,"$dp_dir/$dp_file");
				my $dpcolum=<DP>;chomp $dpcolum;
				my @dpcolum=split(/\t/,$dpcolum);
				for(my $i=2;@dpcolum>$i;$i++){
						if($patient{$dpcolum[$i]}{gender} eq "male"){push(@male,$i);}
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

				
