#!/usr/bin/perl
use strict;
use warnings;

#perl bamslicing.pl body_part

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $bp=$ARGV[0]; #body part
my $bamdir="$bp/$bp"."_norm_cancergenes_bam";
mkdir $bamdir;

my $token_path=`ls ~/git/gdc_il|grep 'gdc-user-token'`;
if($token_path!~/./){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="~/git/gdc_il/$token_path";
my $token=`cat $token_path`;

my $manifest=`ls $bp|grep gdc_manifest`;chomp $manifest;

open(MAN,"$bp/$manifest") or die "manifest file:$manifest cannot open\n";
#my($start_line_num,$end_line_num)=($ARGV[2],$ARGV[3]);
my$line_num=0;
use Parallel::ForkManager;
my $max_processes=10;
my $pm = new Parallel::ForkManager($max_processes);
my $dev_null=<MAN>;
my $json=$ENV{"HOME"}."/git/gdc_il/ilp27/maf_norm/top_driverallexon_json.txt";
($json and -e $json) or die "ERROR:$json is not exist!^n";

while(<MAN>){
		$line_num++;
		if($line_num <345){next;}
		#forks and returns the pid for child
#		my $pid = $pm->start and next;

#		if($start_line_num > $line_num){next;
#		}elsif($end_line_num < $line_num){last;}
		chomp;
		my @line=split(/\t/,);
		print "$line_num:curl $line[1] now\n";
		#system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $bamdir/$line[1] > /dev/null 2>&1");
		system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $bamdir/$line[1]");
		last;

#		$pm->finish; #terminates the child process
}
close MAN;
exit;
