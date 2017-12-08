#!/usr/bin/perl
use strict;
use warnings;

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc"){die "ERROR:doing on wrong dir\nthis perl must doing on /Volumes/areca42TB2/gdc\n";}

my $bp=$ARGV[0]; #body part
my $dpdir="$pwd/varscan/$bp/depth";

#read patient list and gender file
my %gender=();
open(GF,"$pwd/varscan/$bp/gender_age.tsv");
<GF>;
while(<GF>){
		chomp;
		my @line=split(/\t/,);
		$gender{$line[0]}=$line[1];
}
close GF;

my @dp_files=`ls $dpdir|grep out`;chomp @dp_files;
my%ac=();
foreach my $dp_file(@dp_files){
		open(DP,"$dpdir/$dp_file");
		my @dpcolum = split(/\t/,<DP>);
		my @male_cnum=();
		for(my $i=2;@dpcolum>$i;$i++){
				if(!defined $gender{$dpcolum[$i]}){next;
				}elsif($gender{$dpcolum[$i]} eq "male"){
						push(@male_cnum,$i);
				}
		}
		while(<DP>){
				if($_ !~ /^chrX/){next;}
				chomp;
				my @line=split(/\t/,);
				if(!defined$ac{$line[1]}){$ac{$line[1]}=0;}
				foreach my $i (@male_cnum){
						if($line[$i] >= 8){$ac{$line[1]}++;}
				}
		}
		close DP;
}

open(OUT,">$dpdir/coverage_x_male.tsv");
foreach my $posi (sort{ $a<=>$b} keys %ac){
		print OUT "chrX\t$posi\t$ac{$posi}\n";
}
close OUT;
exit;
