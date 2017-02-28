#!/usr/bin/perl
use strict;
use warnings;

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB/tcga"){die "ERROR:doing on wrong dir\nthis perl must doing on /Volumes/areca42TB/tcga\n";}

my $bp=$ARGV[0]; #body part
my $dpdir="$pwd/maf_norm/$bp/depth";

#make body_part sex_json

my ($json,$jsonbp,$response_sexfile,$response_cel,$file_sex,$list_normbam_perfect,$ascat_error)=
	("$pwd/maf_norm/sex_json.txt",
	"$dpdir/$bp"."_sex_json.txt",
	"$dpdir/sex_file.tsv",
	"$pwd/CNA/$bp/cel/response_normcel.tsv",
	"$pwd/CNA/$bp/cel/file_sex",
	"$dpdir/list_of_perfect_maf.tsv",
	"$pwd/CNA/$bp/cel/error_of_annotate_ascat.txt");
#make json file of gender curl
my $bp_for_json=ucfirst $bp;
system("cat $json|sed s/BodyPart/$bp_for_json/ >$jsonbp");

#check need file exist?
if(! -f $json){die "ERROR:$json file is not exist\n";}
if(! -f $response_cel){die "ERROR:$response_cel is not exist\n";}
if(! -f $file_sex){die "ERROR:$file_sex is not exist\n";}
if(! -f $list_normbam_perfect){die "ERROR:$list_normbam_perfect file is not exist\n";}
if(! -f $ascat_error){die "ERROR:$ascat_error file is not exist\n";}

#list up having ascat CNA data
my @list_perf=`ls $pwd/CNA/$bp/cel/ascat/`;chomp @list_perf;
#mix ascat(CNA) & norm_maf patient both ok, so now remove ascat error
open(ER,"$ascat_error");
my $devnull=<ER>;
while(<ER>){
		chomp;
		my $line=$_;
		@list_perf=grep{ $line ne $_ } @list_perf;
}
close ER;

#filter maf_norm too ok
open(LP,"$list_normbam_perfect");
my %lnp=();
$devnull=<LP>;
while(<LP>){
		chomp;
		my @line=split(/\t/,);
		if($line[1] eq "ok"){$lnp{$line[0]}="ok";}
}
close LP;
@list_perf=grep{defined$lnp{$_}} @list_perf;
my %list_perf=map{$_=>"ok"}@list_perf; #@list_perf having patient_id is both ascat CNA and maf_norm ok

#nkf command installed
my $nkfpath=`which nkf`;chomp $nkfpath;
#curl gender response
if($nkfpath and -e $nkfpath){die "ERROR:nkf was not installed. please do\nbrew install nkf\n";}
system("curl --request POST --header \"Content-Type: application/json\" --data \@$jsonbp 'https://gdc-api.nci.nih.gov/cases'|nkf -Lu >$response_sexfile");

#read file_sex(of CNA) and make hash (%cel)
open(CEL,"$response_cel") or die "ERROR:cant open $response_cel\n";
my ($submitter_colum,$file_colum)=("x","x");
my $colum=<CEL>; #first line (= colum name)
chomp $colum;
my @colname=split(/\t/,$colum);
for(my $i=0;@colname>$i;$i++){
		if($colname[$i] eq "cases_0_submitter_id"){$submitter_colum=$i;}
		if($colname[$i] eq "file_name"){$file_colum=$i;}
}
if(($submitter_colum eq "x")||($file_colum eq "x")){die "ERROR:submitter_id or file_name colum is not exist on $response_cel\n";}

my %cel=();
while(<CEL>){
		chomp;
		my @line=split(/\t/,);
		$cel{$line[$file_colum]}=$line[$submitter_colum];
}
close CEL;

#read file_sex made by doing ascat and make patient_id => gender hash(of CNA: %file_sex
my %file_sex=();
open(FS,"$file_sex") or die "ERROR:cant open $file_sex\n";
while(<FS>){
		chomp;
		my @line=split(/\t/,);
		$file_sex{$cel{$line[0]}}=$line[1];
}
close FS;

#read geneder response and out need patient gender_file & hash(%list_perf)
open(SF,"$response_sexfile");
my ($id,$gender)=("x","x");
$colum=<SF>;chomp $colum;
@colname=split(/\t/,$colum);
for(my $i=0;scalar(@colname)>$i;$i++){
		if($colname[$i] eq "submitter_id"){$id=$i;}
		elsif($colname[$i] eq "demographic_gender"){$gender=$i;}
}
if(($id eq "x")||($gender eq "x")){die "ERROR:submitter_id or demographic_gender colum is not exist on $response_sexfile\n";}

open(OUT,">$dpdir/gender_file.tsv");
print OUT "patient_id\tgender\n";

while(<SF>){
		chomp;
		my @line=split(/\t/,);
		if(!defined $list_perf{$line[$id]}){next;}
		if((!defined $line[1])||($line[0] eq "")||($line[$gender] eq "")){# gender no data patient so put array sex data
				print OUT "$line[$id]\t$file_sex{$line[$id]}\n";
				$list_perf{$line[$id]}=$file_sex{$line[$id]};
		}else{
				print OUT"$line[$id]\t$line[$gender]\n";
				$list_perf{$line[$id]}=$line[$gender];
		}
}
close SF;
close OUT;

#make all coverage file
my @dpls=`ls $dpdir|grep out`;chomp @dpls;
my %af=();
foreach my $dpls (@dpls){
		print "reading $dpls\n";
		$|=1;
		my @ok_patient_colum_num=();
		my @X_patient_colum_num=();
		my @Y_patient_colum_num=();
		open(DP,"$dpdir/$dpls");
		my $dpcolum=<DP>;chomp $dpcolum;
		my @dpcolum=split(/\t/,$dpcolum);
		for(my $i=2;@dpcolum>$i;$i++){
				if(!defined $list_perf{$dpcolum[$i]}){next;}
				elsif($list_perf{$dpcolum[$i]} eq "female"){
						push(@ok_patient_colum_num,$i);
						push(@X_patient_colum_num,$i);
				}elsif($list_perf{$dpcolum[$i]} eq "male"){
						push(@ok_patient_colum_num,$i);
						push(@Y_patient_colum_num,$i);
				}else{
						print "WARNING:$dpdir/$dpls colum $i $dpcolum[$i] has no gender data, so we defined the patient is female XX\n";
						push(@ok_patient_colum_num,$i);
						push(@X_patient_colum_num,$i);
				}
		}
		while(<DP>){
				chomp;
				my @line=split(/\t/,);
				$line[0] =~ s/^chr//;
				if(!defined $af{$line[0]}{$line[1]}){$af{$line[0]}{$line[1]}=0;}
				if($line[0] ne "X"){
						foreach my $coln(@ok_patient_colum_num){
								if($line[$coln] > 4){$af{$line[0]}{$line[1]} += 2;}
						}
				}else{
						foreach my $coln(@X_patient_colum_num){
								if($line[$coln] > 4){$af{$line[0]}{$line[1]} += 2;}
						}
						foreach my $coln(@Y_patient_colum_num){
								if($line[$coln] > 4){$af{$line[0]}{$line[1]}++;}
						}
				}
		}
		close DP;
}
#書き出す順番用
my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);
open(OUT,">$dpdir/coverage_all_data_exist_patient.tsv");
foreach my $chr(@chr){
		foreach my $position (sort{$a<=>$b}keys %{$af{$chr}}){
				print OUT "$chr\t$position\t$af{$chr}{$position}\n";
		}
}
close OUT;

exit;

