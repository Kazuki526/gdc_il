#!/usr/vin/perl
use warnings;
use strict;

my $body_part=$ARGV[0];
my @bad=`grep 'redownloaded' $body_part/out_bam2maf.log`;

my @bad_line=();
foreach my $line(@bad){
		chomp $line;
		if($line =~ /^(\d+):.+redownloaded (chr.+):(\d+) is ok$/){
				if(($2 ne "chrX")||($3 < 134428791)){push(@bad_line,$1);}
		}elsif($line =~ /^(\d+):/){push(@bad_line,$1);}
}

open(RES,"/Volumes/cancer/kaz_gdc/$body_part/response.tsv");
my $colum=<RES>;
my @colum=split(/\t/,$colum);
my($patient_id_colum_num,$filename_colum_num)=("","");
for (my $i=0;@colum > $i;$i++){
		if($colum[$i] eq "cases_0_submitter_id"){$patient_id_colum_num=$i;
		}elsif($colum[$i] eq "file_name"){$filename_colum_num=$i;}
}
if(($patient_id_colum_num eq "")||($filename_colum_num eq "")){die "!!!ERRORRRRRR!!! $body_part/response.tsv colum is not normal!\n";}

my%response=();
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if($response{$line[$patient_id_colum_num]}){
				$response{$line[$patient_id_colum_num]}{files}.="\t$line[$filename_colum_num]";
				$response{$line[$patient_id_colum_num]}{filenum}++;
		}else{
				$response{$line[$patient_id_colum_num]}{files} ="$line[$filename_colum_num]";
				$response{$line[$patient_id_colum_num]}{filenum}=1;
		}
}
close RES;

my $manifest=`ls /Volumes/cancer/kaz_gdc/$body_part/|grep 'manifest'`;
chomp $manifest;
open(MANI,"$body_part/$manifest");
my $dev_null=<MANI>;
my %mani_ln=();
my $line_num=0;
while(<MANI>){
		$line_num++;
		chomp;
		my @line=split(/\t/,);
		$mani_ln{$line[1]}=$line_num;
}
close MANI;

open(OUT,">$body_part/list_of_perfect_maf.tsv");
print OUT "patient_id\tfocal\n";
foreach my $patient (keys %response){
		if($response{$patient}{filenum}==1){
				if(grep{$_ == $mani_ln{$response{$patient}{files}}}@bad_line){print OUT "$patient\tno\n";}
				else{print OUT "$patient\tok\n";}
		}else{
				my $focal=0;
				my @files=split(/\t/,$response{$patient}{files});
				foreach my $files (@files){
						if(grep{$_ == $mani_ln{$files}}@bad_line){$focal++;}
				}
				if($focal==$response{$patient}{filenum}){print OUT "$patient\tno\n";}
				else{print OUT "$patient\tok\n";}
		}
}
close OUT;
exit;

