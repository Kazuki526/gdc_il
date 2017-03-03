#!/usr/bin/perl
use warnings;
use strict;

#dir check
my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $token_path=`ls ~/git/gdc_il|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="~/git/gdc_il/$token_path";
my $token=`cat $token_path`;

my $bp=$ARGV[0];
my $bamdir="$bp/$bp"."_norm_cancergenes_bam";
my $manifest = `ls $bp|grep '^gdc_manifest'`;chomp $manifest;

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(10);

open(MANI,"$bp/$manifest");
my $devnull=<MANI>;
my $linen=0;
my $json=$ENV{"HOME"}."/git/gdc_il/ilp27/maf_norm/top_driverallexon_json.txt";
($json and -e $json) or die "ERROR:$json not exist!!";
while(<MANI>){
		$linen++;
		#forks and returns the pid for child
		my $pid = $pm->start and next;
		
		chomp ;
		my @line=split(/\t/,);
		print "$linen:$line[1] doing\n";
		my $taila=`samtools view $bamdir/$line[1] |tail -n 1`;
		if(!$taila){$taila="0\t0\t0";}
		my @taila=split(/\t/,$taila);
		if(($taila[2] ne "chrX")||($taila[3] < 134428791 )){
				print "$linen:$line[1] chr $taila[2]:$taila[3] eq not chrX:134428791~ so download again\n";
				my $focal=0;
				while($focal==0){
						system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@$json --output $bamdir/$line[1] > /dev/null 2>&1");
						my $tailb=`samtools view $bamdir/$line[1] |tail -n 1`;
						my @tailb=split(/\t/,$tailb);
						if(!defined $tailb){$taila="0\t0\t0";}
						if(($tailb[2] ne "chrX")&&($tailb[3] > 134428791)&&($taila ne $tailb)){$taila=$tailb;
						}else{
								$focal++;
								print "$linen:$line[1] redownloaded $tailb[2]:$tailb[3] is ok\n";
								}
				}
		}

		$pm->finish; #terminates the child process
}

