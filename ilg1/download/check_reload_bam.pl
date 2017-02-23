#!/usr/bin/perl
use warnings;
use strict;

#dir check

my $token_file=`ls ~/git/gdc_il/|grep '^gdc-user-token'`;chomp $token_file;
my $token=`cat ~/git/gdc_il/$token_file`;
my $bp=$ARGV[0];
my $manifest = `ls $bp|grep '^gdc_manifest'`;chomp $manifest;

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(10);

open(MANI,"$bp/$manifest");
my $devnull=<MANI>;
my $linen=0;
while(<MANI>){
		$linen++;
		#forks and returns the pid for child
		my $pid = $pm->start and next;
		
		chomp ;
		my @line=split(/\t/,);
		print "$linen:$line[1] doing\n";
		my $taila=`samtools view $bp/$line[1] |tail -n 1`;
		if(!defined $taila){$taila="0\t0\t0";}
		my @taila=split(/\t/,$taila);
		if(($taila[2] ne "chrX")||($taila[3] < 134428791 )){
				print "$linen:$line[1]  $taila[2]:$tailb[3] eq not chrXorY so download again\n";
				my $focal=0;
				while($focal==0){
						system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@top_driverallexon_json.txt --output $bp/$line[1] > /dev/null 2>&1");
						my $tailb=`samtools view $bp/$line[1] |tail -n 1`;
						my @tailb=split(/\t/,$tailb);
						if(!defined $tailb){$taila="0\t0\t0";}
						if(($tailb[2] ne "chrX")&&($tailb[3] > 134428791 )&&($taila ne $tailb)){$taila=$tailb;
						}else{
								$focal++;
								print "$linen:$line[1] redownloaded $tailb[2]:$tailb[3] is ok\n";
								}
				}
		}

		$pm->finish; #terminates the child process
}

