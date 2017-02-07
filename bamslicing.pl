#!/usr/bin/perl
use strict;
use warnings;

#perl bamslicing.pl MANIFEST.tsv out_DIR maxparallel

my $token=`cat gdc-user-token.txt`;
open(MAN,"$ARGV[0]");
#my($start_line_num,$end_line_num)=($ARGV[2],$ARGV[3]);
my$line_num=0;
use Parallel::ForkManager;
my $max_processes=$ARGV[2];
my $pm = new Parallel::ForkManager($max_processes);
my $dev_null=<MAN>;

while(<MAN>){
		$line_num++;
		
		#forks and returns the pid for child
		my $pid = $pm->start and next;

#		if($start_line_num > $line_num){next;
#		}elsif($end_line_num < $line_num){last;}
		chomp;
		my @line=split(/\t/,);
		print "$line_num:curl $line[1] now\n";
		system("curl --header \"X-Auth-Token: $token\" --request POST https://gdc-api.nci.nih.gov/slicing/view/$line[0] --header \"Content-Type: application/json\" -d\@top_driverallexon_json.txt --output $ARGV[1]/$line[1] > /dev/null 2>&1");

		$pm->finish; #terminates the child process
}
exit;
