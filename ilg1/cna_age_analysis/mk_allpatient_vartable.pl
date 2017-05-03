#!/usr/bin/perl
use warnings;
use strict;

my @bp=qw(breast lung brain kidney);

my %pa_inf=();
open(PA,"$ENV{HOME}/git/all_patient/all_patient_response.tsv");
my $colum=<PA>;
while(<PA>){
		chomp;
		my @line=split(/\t/, $_, -1);
		if($line[4] eq ""){next;}
		$pa_inf{$line[1]}{cancer_type}=$line[0];
		$pa_inf{$line[1]}{age}        =$line[4];
		$pa_inf{$line[1]}{die}        =$line[5];
		$pa_inf{$line[1]}{smoke}      =$line[6];
}
close PA;

my %chr_leng=();
open(CHR,"/Volumes/areca42TB/tcga/CNA/all_patient/chr_start_end_pos.tsv");
<CHR>;
while(<CHR>){
		chomp;
		my @line=split(/\t/,);
		$chr_leng{$line[0]}{start}=$line[1];
		$chr_leng{$line[0]}{end}=$line[2];
}
close CHR;

open(OUT,">/Volumes/areca42TB/tcga/CNA/all_patient/all_patient_cna_length.tsv");
print OUT "patient_id\tbody_part\tcancer_type\tdiagnoses_age\tday_to_death\tsmoke_year\tWGD\tamp_chr_num\tleng_0\tleng_1\tleng_2\tleng_3\tleng_4\tleng_5_\tline_num\n";
foreach my $bp (@bp){
		my $ascat_dir="/Volumes/areca42TB/tcga/CNA/$bp/cel/ascat";
		my @ls=`ls $ascat_dir`;chomp @ls;
		open(OUTA,"|gzip -c >/Volumes/areca42TB/tcga/CNA/all_patient/$bp"."_ascat_out.tsv.gz");
		print OUTA "patient_id\tchr\tstartpos\tendpos\tnMajor\tnMinor\tploidy\tpurity\tchr_duplicate\tWG_dupli\n";
		foreach my $pa(@ls){
				if(!defined $pa_inf{$pa}){next;}
				open(IN,"$ascat_dir/$pa/$pa"."_ascat.tsv")or next;
				my $colum=<IN>;
				my @leng=(0,0,0,0,0,0);
				my @wgd=(0,0,0,0,0,0);
				my @chr_d=(0,0,0,0,0,0);
				my @lines=();
				my @all_lines=();
				my ($chr,$wc_line)=("");
				while(<IN>){
						$wc_line++;
						chomp;
						my @line=split(/\t/,);
						if($chr ne $line[1]){
								if($chr eq ""){$chr=$line[1];
								}else{
										my $chr_d=&chr_duplicate($chr,@chr_d);
										$wgd[$chr_d]++;
										for(my $chri=0;5>=$chri;$chri++){$leng[$chri]+=$chr_d[$chri];}
										if($line[1] eq "X"){last;
										}else{
												$chr=$line[1];
												for(my $l=0;@lines>$l;$l++){push(@all_lines,"$lines[$l]\t$chr_d");}
												@lines=();
												@chr_d=(0,0,0,0,0,0);
										}
								}
						}
						push(@lines,$_);
						if(!(($line[4] ==1)&&($line[5]==1))){
								my $n=$line[4] + $line[5];
								$chr_d[$n]+=$line[3] -$line[2] +1;
						}
				}
				close IN;
				my $wgd=&whole_genome_duplicate(@wgd);
				my $amp_chr=$wgd[3] + $wgd[4] + $wgd[5];
				print OUT "$pa\t$bp\t$pa_inf{$pa}{cancer_type}\t$pa_inf{$pa}{age}\t$pa_inf{$pa}{die}\t$pa_inf{$pa}{smoke}\t$wgd\t$amp_chr\t".join("\t",@leng)."\t$wc_line\n";
				print OUTA join("\t$wgd\n",@all_lines)."\t$wgd\n";
		}
		close OUTA;
}
close OUT;



#=========================================================================================================
sub chr_duplicate {
		my ($chr,@chr_dupli)=@_;
		my $chr_leng=$chr_leng{$chr}{end} - $chr_leng{$chr}{start};
		if($chr_dupli[4]+$chr_dupli[5] >$chr_leng *0.6){
				return(4);
		}elsif($chr_dupli[3] + $chr_dupli[4] + $chr_dupli[5] >$chr_leng*0.6){
				return(3);
		}else{return(2);}
}
#========================================================================================================
sub whole_genome_duplicate{
		my @wgd_=@_;
		if($wgd_[4] + $wgd_[5] > $wgd_[3] + $wgd_[2]){
				return(4);
		}elsif($wgd_[3] + $wgd_[4] +$wgd_[5] > $wgd_[2]){
				return(3);
		}else{return(2);}
}

