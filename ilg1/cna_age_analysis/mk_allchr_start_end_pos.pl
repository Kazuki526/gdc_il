#!/usr/bin/perl
use strict;
use warnings;


#brain lungではchr11のstartが188572(他は188510)だが、ほぼ誤差なのでbreastのデータを使用することにした
#これを実行後 
#mv /Volumes/areca42TB/tcga/CNA/all_patient/breast_chr.tsv /Volumes/areca42TB/tcga/CNA/all_patient/chr_start_end_pos.tsv
my $bp="breast";
my $ascat_dir="/Volumes/areca42TB/tcga/CNA/$bp/cel/ascat";
my @ls=`ls $ascat_dir`;chomp @ls;

my @start=();
my @end=();
for(my $i=0;100>$i;$i++){
		my $pa=$ls[$i];
		open(IN,"$ascat_dir/$pa/$pa"."_ascat.tsv") or next;
		my $colum=<IN>;
		my %chr=();
		while(<IN>){
				chomp;
				my @line=split(/\t/,);
				if($line[1] eq "X"){last;}
				if(!defined$chr{$line[1]}){
						$chr{$line[1]}{start}=$line[2];
						$chr{$line[1]}{end}=$line[3];
				}else{
						if($line[2] < $chr{$line[1]}{start}){$chr{$line[1]}{start}=$line[2];}
						if($line[3] > $chr{$line[1]}{end})  {$chr{$line[1]}{end}  =$line[3];}
				}
		}
		close IN;
		for(my $chr=1;23>$chr;$chr++){
				if((!defined $chr{$chr}{start})||(!defined$chr{$chr}{end})){
						print "$pa has not have chr$chr\n";
				}
				if($i==0){$start[$chr]=$chr{$chr}{start};$end[$chr]=$chr{$chr}{end};
				}else{
						if($start[$chr] > $chr{$chr}{start}){
								print "$pa has forward start chr$chr\n";
								$start[$chr]=$chr{$chr}{start};
						}
						if($end[$chr] < $chr{$chr}{end}){
								print "$pa has posterior end chr$chr\n";
								$end[$chr]=$chr{$chr}{end};
						}
				}
		}
}

open(OUT,">/Volumes/areca42TB/tcga/CNA/all_patient/$bp"."_chr.tsv");
print OUT "chr\tstart\tend\n";
for(my $chr=1;23>$chr;$chr++){
		print OUT "$chr\t$start[$chr]\t$end[$chr]\n";
}
close OUT;

exit;
