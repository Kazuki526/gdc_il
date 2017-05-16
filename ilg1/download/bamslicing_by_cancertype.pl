#!/usr/bin/perl
use strict;
use warnings;

#dir chech
my $pwd=`pwd`;chomp $pwd;
unless(($pwd eq "/Volumes/areca42TB2/gdc/tumor_bam")||($pwd eq "/Volumes/areca42TB2/gdc/norm_bam")) {die "ERROR: doing on wrong dir !!\n";}

#perl bamslicing.pl MANIFEST.tsv out_DIR maxparallel
my $token_path=`ls ~/git/gdc_il|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="~/git/gdc_il/$token_path";
my $token=`cat $token_path`;

my $ct=$ARGV[0]; #cancer type
my $bamdir=$ct;
my $manifest=`ls $bamdir|grep '^gdc_manifest'`;chomp $manifest;
($manifest and -e "$bamdir/$manifest") or die "ERROR: there arenot manifest file!!\n";

open(MAN,"$bamdir/$manifest");
#my($start_line_num,$end_line_num)=($ARGV[2],$ARGV[3]);
my $linen=0;
use Parallel::ForkManager;
my $max_processes=10;
my $pm = new Parallel::ForkManager($max_processes);
my $dev_null=<MAN>;
my $json=$ENV{"HOME"}."/git/driver_genes/onlytop105/top_driver105exon_json.txt";
(-e $json) or die "ERROR:$json not exist!!\n";
open(ERROR,">$bamdir/download_errored_file.txt");

while(<MAN>){
		$linen++;
		
		#forks and returns the pid for child
		my $pid = $pm->start and next;

#		if($start_line_num > $line_num){next;
#		}elsif($end_line_num < $line_num){last;}
		chomp;
		my @line=split(/\t/,);
		my @ls=`ls $bamdir`;chomp @ls;
		if(!grep{$_ eq "$line[1]"}@ls){
				print "$linen:curl $line[1] now\n";
				system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $bamdir/$line[1] > /dev/null 2>&1");
		}
		print "$linen:$line[1] checking\n";
		my $taila;
		my $focal=0;
		my $time=0;
		while($focal==0){
				$time++;
				my $tailb=`samtools view $bamdir/$line[1] |tail -n 1`;
				if(!$tailb){$tailb="0\t0\t0\t0";}
				my @tailb=split(/\t/,$tailb);
				if($time >10){
						print "$linen:$line[1] download more than 10 times so this file cannot download? redownload result $tailb[2]:$tailb[3]\n";
						$focal++;
						print ERROR "$line[1]\n";
				}elsif(($tailb[2] eq "chrX")&&($tailb[3] > 134428789)){
						$focal++;
						print "$linen:$line[1] redownloaded $tailb[2]:$tailb[3] is ok\n";
				}elsif(($tailb[3] != 0)&&(`samtools view $bamdir/$line[1] 2>&1|head -n 1` !~ /EOF\smarker\sis\sabsent/)){
						$focal++;
						print "$linen:$line[1] redownloaded $tailb[2]:$tailb[3] is ok\n";
				}else{
						$taila=$tailb;
						print "$linen:$line[1] download error. download again\n";
				}
				if($focal==0){
						system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $bamdir/$line[1] > /dev/null 2>&1");
				}
		}

		$pm->finish; #terminates the child process
}
close ERROR;
exit;
