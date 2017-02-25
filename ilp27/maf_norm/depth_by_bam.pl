#!/usr/bin/perl
use warnings;
use strict;

#perl depth_by_bam.pl body_part 

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $bp=$ARGV[0]; #body part
my $bamdir="$bp/$bp"."_norm_cancergenes_bam";

system("perl ~/git/gdc_il/ilp27/maf_norm/make_list_of_perfect_maf.pl $bp");

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
#				print "$response{$pid}\n";
				$response{$pid}="$pid.bam";
		}else{
				$response{$pid} =~ s/;//;
#				print"1bam $file[0]\n";
		}
}
mkdir "$bp/depth" or print "depth dir is already exist!!\n";

open(LIST,"./$bp/list_of_perfect_maf.tsv");
<LIST>;
my %list=();
while(<LIST>){
		chomp;
		my @line=split(/\t/,);
		$list{$line[0]}=$line[1];
}
close LIST;

my $file_num=0;
my @bam=();
my @pid=();
foreach my $pid(keys %response){
		if($list{$pid} eq "no"){next;}
		push(@pid,$pid);
		push(@bam,"$bamdir/$response{$pid}");
		if(scalar(@bam) ==50){
				$file_num++;
				print "$file_num th file writing\n";
				my $outfile="out$file_num.tsv";
				open(OUT,">$bp/depth/$outfile");
				print OUT "chr\tposition\t".join("\t",@pid)."\n";
				close OUT;
				`samtools depth -q 13 -b top_driver105.bed @bam >>$bp/depth/$outfile`;
				@bam=();
				@pid=();
		}
		
}
$file_num++;
print "$file_num th file writing\n";
my $outfile="out$file_num.tsv";
open(OUT,">$bp/depth/$outfile");
print OUT "chr\tposition\t".join("\t",@pid)."\n";
close OUT;
`samtools depth -q 13 -b top_driver105.bed @bam >>$bp/depth/$outfile`;

exit;
