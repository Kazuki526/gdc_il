#!/usr/bin/perl
use warnings;
use strict;

#doing on /Volumes/areca42TB/gdc
#perl ~/git/gdc_il/varscan/varscan_genderage.pl bodypart 2>&1|tee varscan/bodypart/out.log


my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc"){die "ERROR:doing on wrong dir!\n cd /Volumes/areca42TB2/gdc\n";}

#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

my $pj=$ARGV[0];#project_id

#check reference of bam are exist?
my $ref="GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";

#check top driver 105genes bed file exist?
my $bed="/Users/kaz/git/gdc_il/ilp27/maf_norm/top_driver105.bed";
(-e $bed) or die "ERROR:not exist bed file:$bed\n";

#check pj ascat file exist ?
my $ascat_file="/Volumes/areca42TB/tcga/CNA/$pj/cel/annotate_ascat.tsv.gz";
(-e $ascat_file)or die "ERROR:not exist ascat file:$ascat_file\n";

#check download errored file
my @tumor_error=();
my @norm_error=();
open(TE,"tumor_bam/$pj/download_errored_file.txt") or die "ERROR:there is not exist tumor download error list file download_errored_file.txt!!\n";
@tumor_error=<TE>;chomp @tumor_error;
close TE;
open(NE,"tumor_bam/$pj/download_errored_file.txt") or die "ERROR:there is not exist normal download error list file download_errored_file.txt!!\n";
@norm_error=<NE>;chomp @norm_error;
close NE;

my ($json_for_norm,$json_for_tumor)=("normbam_prjson.txt","tumorbam_prjson.txt");
(-e $json_for_norm and -e $json_for_tumor) or die "ERROR:json files are not exist $pwd/$json_for_norm or $pwd/$json_for_tumor\n";
my $pj_json= uc $pj;
my($norm_json,$tumor_json)=("varscan/$pj/norm.json","varscan/$pj/tumor.json");
`cat $json_for_norm|sed s/cancertype/$pj_json/ >$norm_json`;
`cat $json_for_tumor|sed s/cancertype/$pj_json/ >$tumor_json`;
my($norm_response,$tumor_response)=("varscan/$pj/response_norm.tsv","varscan/$pj/response_tumor.tsv");
print "curl normal and tumor response files\n";
`curl --request POST --header \"Content-Type: application/json\" --data \@$norm_json 'https://gdc-api.nci.nih.gov/files'|nkf -Lu >$norm_response`;
`curl --request POST --header \"Content-Type: application/json\" --data \@$tumor_json 'https://gdc-api.nci.nih.gov/files'|nkf -Lu >$tumor_response`;

my %patient_inf;
print "reading norm response\n";
open(RN,"$norm_response") or die "ERROR:cant open norm response file!!\n";
my $colums=<RN>;chomp $colums;
my @colum=split(/\t/,$colums);
my($file_name,$gender,$age,$case_id)=("","","","");
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "file_name"){$file_name=$i;
		}elsif($colum[$i] eq "cases_0_demographic_gender"){$gender=$i;
		}elsif($colum[$i] eq "cases_0_submitter_id"){$case_id=$i;
		}elsif($colum[$i] eq "cases_0_diagnoses_0_age_at_diagnosis"){$age=$i;
		}
}
($file_name ne "" and $gender ne "" and $case_id ne "" and $age ne "") or die "ERROR:response norm file's colum is not normaly:$file_name,$gender,$case_id,$age\n";

while(<RN>){
		chomp;
		my @line=split(/\t/,);
		if(scalar @line != 5){
				for(my$i=0;5>=$i;$i++){
						if(!defined$line[$i] ){$line[$i] ="";}
				}
		}
		if(defined $patient_inf{$line[$case_id]}){
				$patient_inf{$line[$case_id]}{'file_norm'}.=";$line[$file_name]";
				if($patient_inf{$line[$case_id]}{'gender'} ne $line[$gender] and $patient_inf{$line[$case_id]}{'gender'} and $line[$gender]){
						print "WARNING:$line[$case_id] norm have diffrent gender\n";}
				if($patient_inf{$line[$case_id]}{'age'} ne $line[$age] and $patient_inf{$line[$case_id]}{'age'} and $line[$age]){
						print "WARNING:$line[$case_id] norm have diffrent age\n";}
		}else{
				$patient_inf{$line[$case_id]}{'file_norm'}=$line[$file_name];
				$patient_inf{$line[$case_id]}{'age'}      =$line[$age];
				$patient_inf{$line[$case_id]}{'gender'}   =$line[$gender];
		}
}
close RN;

print "reading tumor response\n";
open(RT,"$tumor_response") or die "ERROR:cant open tumor response file!!\n";
$colums=<RT>;chomp $colums;
@colum=split(/\t/,$colums);
($file_name,$gender,$age,$case_id)=("","","","");
for(my $i=0;@colum>$i;$i++){
		if($colum[$i] eq "file_name"){$file_name=$i;
		}elsif($colum[$i] eq "cases_0_demographic_gender"){$gender=$i;
		}elsif($colum[$i] eq "cases_0_submitter_id"){$case_id=$i;
		}elsif($colum[$i] eq "cases_0_diagnoses_0_age_at_diagnosis"){$age=$i;
		}
}
($file_name ne "" and $gender ne "" and $case_id ne "" and $age ne "") or die "ERROR:response tumor file's colum is not normaly\n";

while(<RT>){
		chomp;
		my @line=split(/\t/,);
		if(scalar @line != 5){
				for(my$i=0;5>=$i;$i++){
						if(!defined$line[$i] ){$line[$i] ="";}
				}
		}
		if(defined $patient_inf{$line[$case_id]}){
				if(defined $patient_inf{$line[$case_id]}{'file_tumor'}){
						$patient_inf{$line[$case_id]}{'file_tumor'}.=";$line[$file_name]";
				}else{
						$patient_inf{$line[$case_id]}{'file_tumor'}=$line[$file_name];
				}
				if($patient_inf{$line[$case_id]}{'gender'} ne $line[$gender] and $patient_inf{$line[$case_id]}{'gender'} and $line[$gender]){
						print "WARNING:$line[$case_id] have diffrent gender\n";}
				if($patient_inf{$line[$case_id]}{'age'} ne $line[$age] and $patient_inf{$line[$case_id]}{'age'} and $line[$age]){
						print "WARNING:$line[$case_id] have diffrent age\n";}
		}
}
close RT;

print "reading ascat files\n";
my %ascat;
open(AF,"gunzip -c $ascat_file|")or die "ERROR:cannot open $ascat_file\n";
my$devnull=<AF>;#0:gene_symbol 1:patient_id 2:nmajor 3:nminor 4:chr 5:cna_start 6:cna_end 7:start_rate  8:end_rate 9:gene_start 10:gene_end 11:ploidy 12:purity
my $patient_id="";
while(<AF>){
		chomp;
		my @line=split(/\t/,);
		if($line[1] eq $patient_id){next;}
		$patient_id=$line[1];
		$ascat{$patient_id}{'purity'}=$line[12];
}
close AF;

#check both norm&tumor bam exist or not
my @ok_patient=();
foreach my $pid(keys %patient_inf){
		(defined $patient_inf{$pid}{'file_norm'} and defined $patient_inf{$pid}{'file_tumor'} and defined $ascat{$pid}) or next;
		if(grep{$_ eq $patient_inf{$pid}{file_norm}}@norm_error){next;}
		if(grep{$_ eq $patient_inf{$pid}{file_tumor}}@tumor_error){next;}
		if($patient_inf{$pid}{'file_norm'} =~ /;/){
				chdir "norm_bam/$pj";
				my @bam=split(/;/,$patient_inf{$pid}{'file_norm'});
				`samtools merge -f $pid.bam @bam`;
				`samtools index $pid.bam`;
				chdir $pwd;
				$patient_inf{$pid}{'file_norm'}="$pid.bam";}
		if($patient_inf{$pid}{'file_tumor'} =~ /;/){
				chdir "tumor_bam/$pj";
				my @bam=split(/;/,$patient_inf{$pid}{'file_tumor'});
				`samtools merge -f $pid.bam @bam`;
				`samtools index $pid.bam`;
				chdir $pwd;
				$patient_inf{$pid}{'file_tumor'}="$pid.bam";
		}
		push(@ok_patient,$pid);
}

mkdir "varscan/$pj/out";
open(OUT,">varscan/$pj/gender_age.tsv");
print OUT "patient_id\tgender\tage\n";
foreach my $pid(@ok_patient){#全patientでvarscan終わる前にgender age だけ先に書き出す
		print OUT "$pid\t$patient_inf{$pid}{gender}\t$patient_inf{$pid}{age}\n";
}
close OUT;
foreach my $pid(@ok_patient){
		print "doing varsca $pid\n";
		my $normpile = "samtools mpileup -q 10 -f $ref norm_bam/$pj/$patient_inf{$pid}{file_norm}";
		my $tumorpile= "samtools mpileup -q 10 -f $ref tumor_bam/$pj/$patient_inf{$pid}{file_tumor}";
		`zsh -c \"varscan somatic <\($normpile\) <\($tumorpile\) varscan/$pj/out/$pid --tumor-purity $ascat{$pid}{purity} --p-value 0.1 --output-vcf 1\"`;
}
exit;

