#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);

my @info_ac=qw(AC AC_Adj AC_Het AC_Hom);
my @info_an=qw(AN AN_Adj);
foreach my $inf(qw(AC Het Hom)){
				foreach my $race(qw(AFR AMR EAS FIN NFE OTH SAS)){
										push(@info_ac,"$inf"."_$race");
												}
}

map{push(@info_an,"AN_".$_)}qw(AFR AMR EAS FIN NFE OTH SAS);
my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);

open (OUTV,">/Volumes/cancer/exac/exac_nontcga_liftovered.vcf");
open (OUTT,"|gzip -c >/Volumes/cancer/exac/exac/nontcga_liftovered_checked_likevcf.tsv.gz");
print OUTT "chr\tstart\tref\talt\t".join("\t",(@info_ac,@info_an))."\n";
foreach my $chr (@chr){
		print "liftover $chr\n";
		my %remap = map{chomp;my @c=split("\t");$c[0]=~s/^chr//;($c[3],"$c[0]:$c[1]-$c[2]")}`~/vep/samtools/bin/liftOver /Volumes/cancer/exac/before_liftover_bed/$chr.bed ~/.vep/hg19ToHg38.over.chain /dev/stdout /dev/null 2>/dev/null`;
		print "printing out $chr \n";
		open (IN,"gunzip -c /Volumes/cancer/exac/before_liftover/$chr.vcf.gz|");
		while(<IN>){
				if($_=~/^#/){
						if($chr eq "1"){print OUT "$_";}
						next;
				}else{
						chomp;
						my @line=split(/\t/,);
						my %INFO=map{my@inf=split(/=/,);($inf[0],$inf[1])}split(/;/,$line[7]);
						my ($chr,$pos,$end_pos,$refseq);
						my $end=$line[1]+1;
						if($remap{"$line[0]:$line[1]-$end"} and $remap{"$line[0]:$line[1]-$end"}=~ /^([^:]+):(\d+)-(\d+)$/){
								($chr,$pos,$end_pos)=($1,$2,$3);
						}else{print  "WARNIG:there cant liftover $_\tso we dont use this variation\n";next;
						}
						my $fasta=`samtools faidx /Users/kazuki/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz $chr:$pos-$pos`;
						if($fasta =~ /\n([ATGC])\n$/){$refseq=$1;
						}else{die "ERROR: what happen?\n$_\n$fasta";}
						if($refseq eq $line[3]){
								if($line[4] =~ /,/){
										my @alt=split(/,/,$line[4]);
										for(my $i=0;@alt>$i;$i++){
												my @info_out=();my @info_tsv=();
												map{my@inf_ac=split(/,/,$INFO{$_});push(@info_out,"$_=$inf_ac[$i]");push(@info_tsv,$inf_ac[$i])}@info_ac;
												map{push(@info_out,"$_=$INFO{$_}");push(@info_tsv,$INFO[$i])}@info_an;
												print OUTV "$chr\t$pos\t$line[2]\t$line[3]\t$alt[$i]\t$line[5]\t$line[6]\t".join(";",@info_out)."\n";
												print OUTT "$chr\t$pos\t$line[3]\t$alt[$i]\t".join("\t",@info_tsv)."\n";
										}
								}else{
										print OUTV "$chr\t$pos\t".join("\t",@line[2..7])."\n";
										my $info_out=join("\t",map{"$INFO{$_}"}(@info_ac,@info_an));
										print OUTT "$chr\t$pos\t$line[3]\t$line[4]\t$info_out\n";
								}
						}else{
								if($line[4] =~ /,/){
										my @alt=split(/,/,$line[4]);
										for(my $i=0;@alt>$i;$i++){
												if($alt[$i] ne $refseq){
														my @info_out=();my @info_tsv=();
														map{my@inf_ac=split(/,/,$INFO{$_});push(@info_out,"$_=$inf_ac[$i]");push(@info_tsv,$inf_ac[$i])}@info_ac;
														map{push(@info_out,"$_=$INFO{$_}");push(@info_tsv,$INFO[$i])}@info_an;
														print OUTV "$chr\t$pos\t$line[2]\t$refseq\t$alt[$i]\t$line[5]\t$line[6]\t".join(";",@info_out)."\n";
														print OUTT "$chr\t$pos\t$refseq\t$alt[$i]\t".join("\t",@info_tsv)."\n";
												}
										}
										my @info_tsv=();
										foreach my $race(("","_Adj","_AFR", "_AMR", "_EAS", "_FIN", "_NFE", "_OTH", "_SAS")){
												my $ac=sum(split(/,/,$INFO{"AC$race"}));
												my $an=$INFO{"AN$race"};
												push(@info_tsv,$an-$ac);
												if($race ne ""){
														push(@info_tsv,("NA","NA"));
												}
										}
										map{push(@info_tsv,$INFO[$i])}@info_an;
										my $t=0;
										my @info_out=();
										foreach my $inf((@info_ac,@info_an)){
												push(@info_out,"$inf=$info_tsv[$t]");
												$t++;
										}
										print OUTV "$chr\t$pos\t$line[2]\t$refseq\t$line[3]\t$line[5]\t$line[6]\t".join(";",@info_out)."\n";
										print OUTT "$chr\t$pos\t$refseq\t$line[3]\t".join("\t",@info_tsv)."\n";
								}else{
										if($line[4] ne $refseq){print "liftover correct?at $_\nliftover:".$remap{"$line[0]:$line[1]"}."\n";}
										my @info_tsv=();
										foreach my $race(("","_Adj","_AFR", "_AMR", "_EAS", "_FIN", "_NFE", "_OTH", "_SAS")){
												my $ac=sum(split(/,/,$INFO{"AC$race"}));
												my $an=$INFO{"AN$race"};
												push(@info_tsv,$an-$ac);
												if($race ne ""){
														push(@info_tsv,("NA","NA"));
												}
										}
										map{push(@info_tsv,$INFO[$i])}@info_an;
										my $t=0;
										my @info_out=();
										foreach my $inf((@info_ac,@info_an)){
												push(@info_out,"$inf=$info_tsv[$t]");
												$t++;
										}
										print OUTV "$chr\t$pos\t$line[2]\t$refseq\t$line[3]\t$line[5]\t$line[6]\t".join(";",@info_out)."\n";
										print OUTT "$chr\t$pos\t$refseq\t$line[3]\t".join("\t",@info_tsv)."\n";
								}
						}
				}
		}
		close IN;
}
close OUTV;close OUTT




open (OUTV,">/Volumes/cancer/exac/exac_nontcga_liftovered_indel.vcf");
open (ERR,"|gzip -c >/Volumes/cancer/exac/exac_nontcga_liftovererror_indel.vcf.gz");
open (OUTT,"|gzip -c >/Volumes/cancer/exac/exac/nontcga_liftovered_checked_likevcf_indel.tsv.gz");
print OUTT "chr\tstart\tref\talt\t".join("\t",(@info_ac,@info_an))."\n";
foreach my $chr (@chr){
		print "liftover $chr\n";
		my %remap = map{chomp;my @c=split("\t");$c[0]=~s/^chr//;($c[3],"$c[0]:$c[1]-$c[2]")}`~/vep/samtools/bin/liftOver /Volumes/cancer/exac/before_liftover_bed_indel/$chr.bed ~/.vep/hg19ToHg38.over.chain /dev/stdout /dev/null 2>/dev/null`;
		print "printing out $chr \n";
		open (IN,"gunzip -c /Volumes/cancer/exac/before_liftover_indel/$chr.vcf.gz|");
		while(<IN>){
				if($_=~/^#/){
						if($chr eq "1"){print OUT "$_";}
						next;
				}else{
						chomp;
						my @line=split(/\t/,);
						my %INFO=map{my@inf=split(/=/,);($inf[0],$inf[1])}split(/;/,$line[7]);
						my ($chr,$pos,$end_pos,$refseq);
						my $end=$line[1]+length($line[3]);
						if($remap{"$line[0]:$line[1]-$end"} and $remap{"$line[0]:$line[1]-$end"}=~ /^([^:]+):(\d+)-(\d+)$/){
								($chr,$pos,$end_pos)=($1,$2,$3 - 1);
						}
						my $fasta=`samtools faidx /Users/kazuki/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz $chr:$pos-$end_pos`;
						my ($region,$refseq)=split(/\n/,$fasta);
						if($refseq eq $line[3]){
								print OUTV "$chr\t$pos\t".join("\t",@line[2..7])."\n";
								my($ref,$alt,$ref_length,$alt_length)=($line[3],$line[4],length($line[3]),length($line[4]));
								while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
										($ref,$alt) = map{$_ =substr($_,1);($_?$_:"-")}($ref,$alt);
										--$ref_length;--$alt_length;++$pos;
								}
								my $info_out=join("\t",map{"$INFO{$_}"}(@info_ac,@info_an));
								print OUTT "$chr\t$pos\t$ref\t$alt\t$info_out\n";
						}else{
								print ERR "$_\n";
						}
				}
		}
		close IN;
}





exit;
