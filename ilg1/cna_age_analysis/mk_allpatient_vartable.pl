#!/usr/bin/perl
use warnings;
use strict;

my @bp=qw(breast lung brain kidney colorectal);

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

open(OUT,">/Volumes/areca42TB/tcga/CNA/all_patient/all_patient_cna_length2.tsv");
print OUT "patient_id\tbody_part\tcancer_type\tdiagnoses_age\tday_to_death\tsmoke_year\tWGD_major\tWGD_minor\tmajor_d3\tmajor_d2\tmajor_d1\tmajor_a1\tmajor_a2\tmajor_a3\tmajor_a4_\tminor_d3\tminor_d2\tminor_d1\tminor_a1\tminor_a2\tminor_a3\tminor_a4_\tline_num\n";
foreach my $bp (@bp){
		my $ascat_dir="/Volumes/areca42TB/tcga/CNA/$bp/cel/ascat";
		my @ls=`ls $ascat_dir`;chomp @ls;
		open(OUTA,"|gzip -c >/Volumes/areca42TB/tcga/CNA/all_patient/$bp"."_ascat_out.tsv.gz");
		print OUTA "patient_id\tchr\tstartpos\tendpos\tnMajor\tnMinor\tploidy\tpurity\tchr_dupli_major\tchr_dupli_minor\tWGD_major\tWGD_minor\n";
		foreach my $pa(@ls){
				if(!defined $pa_inf{$pa}){next;}
				open(IN,"$ascat_dir/$pa/$pa"."_ascat.tsv")or next;
				my $colum=<IN>;
				my @leng_major=(0,0,0,0,0,0,0,0);
				my @leng_minor=(0,0,0,0,0,0,0,0);
				my @wgd_major=(0,0,0,0,0);
				my @wgd_minor=(0,0,0,0,0);
				my @chr_d_major=(0,0,0,0,0);
				my @chr_d_minor=(0,0,0,0,0);
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
#whole chromosome is duplicated or not?
										my $chr_d_major=&chr_duplicate($chr,@chr_d_major);
										my $chr_d_minor=&chr_duplicate($chr,@chr_d_minor);
										$wgd_major[$chr_d_major]++;$wgd_minor[$chr_d_minor]++;
										if($line[1] eq "X"){last;
										}else{
												$chr=$line[1];
												for(my $l=0;@lines>$l;$l++){push(@all_lines,"$lines[$l]\t$chr_d_major\t$chr_d_minor");}
												@lines=();
												@chr_d_major=(0,0,0,0,0);
												@chr_d_minor=(0,0,0,0,0);
										}
								}
						}
						push(@lines,$_);
						my($nmajor,$nminor)=($line[4] >6? 7: $line[4], $line[5] >6? 7: $line[5]);
						$leng_major[$nmajor] += $line[3] - $line[2] +1;
						$leng_minor[$nminor] += $line[3] - $line[2] +1;
						($nmajor,$nminor)=($nmajor >3? 4: $nmajor, $nminor >3? 4: $nminor);
						$chr_d_major[$nmajor] += $line[3] - $line[2] +1;
						$chr_d_minor[$nminor] += $line[3] - $line[2] +1;
				}
				close IN;
				my $wgd_major = &whole_genome_duplicate(@wgd_major);
				my $wgd_minor = &whole_genome_duplicate(@wgd_minor);
				my @out_major = &ampdel_length($wgd_major,@leng_major);
				my @out_minor = &ampdel_length($wgd_minor,@leng_minor);
				print OUT "$pa\t$bp\t$pa_inf{$pa}{cancer_type}\t$pa_inf{$pa}{age}\t$pa_inf{$pa}{die}\t$pa_inf{$pa}{smoke}\t$wgd_major\t$wgd_minor\t".join("\t",@out_major)."\t".join("\t",@out_minor)."\t$wc_line\n";
				print OUTA join("\t$wgd_major\t$wgd_minor\n",@all_lines)."\t$wgd_major\t$wgd_minor\n";
		}
		close OUTA;
}
close OUT;



#=========================================================================================================
sub chr_duplicate {
		my ($chr,@chr_dupli)=@_;
		my $chr_leng=$chr_leng{$chr}{end} - $chr_leng{$chr}{start};
		if($chr_dupli[4] >$chr_leng *0.6){
				return(4);
		}elsif($chr_dupli[3] + $chr_dupli[4] >$chr_leng*0.6){
				return(3);
		}elsif($chr_dupli[2] + $chr_dupli[3] + $chr_dupli[4] >$chr_leng*0.6){
				return(2);
		}else{return(1);
		}
}
#========================================================================================================
sub whole_genome_duplicate{
		my @wgd=@_;
		if($wgd[4] > $wgd[3] + $wgd[2] + $wgd[1] + $wgd[0] ){
				return(4);
		}elsif($wgd[4] + $wgd[3] > $wgd[2] + $wgd[1] + $wgd[0] ){
				return(3);
		}elsif($wgd[4] + $wgd[3] + $wgd[2] > $wgd[1] + $wgd[0] ){
				return(2);
		}else{return(1);
		}
}
#========================================================================================================
sub ampdel_length{
		my($wgd,@leng)=@_;
		my @amp_del_leng=(0,0,0,0,0,0,0,0); #0;del3, 1:del2, 2:del1, 3:amp1, 4:amp2,,,6:amp4~
		for(my $nchr=0;7>=$nchr;$nchr++){
				my $amp_del_posi = $nchr - $wgd +3 >7? 7: $nchr - $wgd +3;
				$amp_del_posi = $amp_del_posi <0? 0: $amp_del_posi;
				$amp_del_leng[$amp_del_posi] += $leng[$nchr];
		}
		splice(@amp_del_leng, 3, 1);
		return(@amp_del_leng);
}




