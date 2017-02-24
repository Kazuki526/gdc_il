#!/usr/bin/perl
use warnings;
use strict;

#perl ~/git/gdc_il/ilp27/maf_norm/bam2maf.pl body_part 2>&1|tee body_part/out_bam2maf.log

#use Parallel::ForkManager;
#my $pm = new Parallel::ForkManager(2);
my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $bp=$ARGV[0]; #body part
my $bamdir="$bp/$bp"."_norm_cancergenes_bam";

system("perl ~/git/gdc_il/ilp27/maf_norm/check_reload_bam.pl $bp");
$|=1;

my $json=`ls $bp|grep 'bam_json.txt'`;chomp $json;
my $json_path="$bp/$json";
print "download response $json_path file\n";
system("curl --request POST --header \"Content-Type: application/json\" --data \@$json_path 'https://gdc-api.nci.nih.gov/files'|nkf -Lu >$bp/response.tsv");

print "reading response file\n";
open(RES,"$bp/response.tsv") or die "!!!ERRORRRRR!!!  $bp/response.tsv isnot exist!!\n";
my $colum=<RES>;
my @colum=split(/\t/,$colum);
my($patient_id_colum_num,$filename_colum_num)=("","");
for (my $i=0;@colum > $i;$i++){
		if($colum[$i] eq "cases_0_submitter_id"){$patient_id_colum_num=$i;
		}elsif($colum[$i] eq "file_name"){$filename_colum_num=$i;}
}
if(($patient_id_colum_num eq "")||($filename_colum_num eq "")){die "!!!ERRORRRRRR!!! $bp/response.tsv colum is not normal!\n";}

my%response=();
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if($response{$line[$patient_id_colum_num]}){
				$response{$line[$patient_id_colum_num]}.="$line[$filename_colum_num];";
		}else{
				$response{$line[$patient_id_colum_num]} ="$line[$filename_colum_num];";
		}
}
close RES;

#merge bam which same patient have more than 2
foreach my $pid ( keys %response ){
		my @file=split(/;/,$response{$pid});
		if(@file > 1 ){
				print "bam merge $pid\n";
				for(my $i=0;@file > $i;$i++){
						$file[$i]="$bamdir/$file[$i]";
				}
				system("samtools merge $bamdir/$pid.bam @file");
				system("samtools index $bamdir/$pid.bam");
				$response{$pid}="$pid.bam";
		}else{
				print "index $pid\n";
				system("samtools index $bamdir/$file[0]");
				$response{$pid} =~ s/;//;
#				print"1bam $file[0]\n";
		}
}
mkdir "$bp/samtools";
mkdir "$bp/maf";

my $line_num=0;
foreach my $pid(keys %response){
		$line_num++;
#forks nad returns the pid for child
#		my $pid = $pm ->start and next;
		
		print "$line_num\tdoing $pid\n";
		system("zsh ~/git/gdc_il/ilp27/maf_norm/top105cancergene_maf.sh $bp $response{$pid} $pid 2>&1 >/dev/null");

#		$pm->finish;
#fork finish
}

exit;
