#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);


my $exac_path=$ENV{"HOME"}."/.vep/ExAC_nonTCGA.r1.sites.vep.vcf.gz";
(-e $exac_path) or die "ERROR:$exac_path is not exist!\n";
my $driver_genes_path="/Volumes/cancer/kaz_gdc/driver_genes.tsv";
(-e $driver_genes_path) or die "ERROR:$driver_genes_path is not exist!\n";

#list up driver gene to %driver_genes
open(DG,"$driver_genes_path");
<DG>; #delete colum line
my %driver_genes=();
while(<DG>){
		chomp;
		my @line=split(/\t/,);
		if($line[1] <4){last;}
		$driver_genes{$line[0]}="ok";
}
close DG;

open(VCF,"gunzip -c $exac_path|");
my %info=();
my @chr;
my ($chr,$comment)=("","");
while(<VCF>){
		if($_ =~ /^##INFO=<ID=CSQ,.+Format:\s(.+)">$/){
				my @info=split(/\|/,$1);
				for(my $i=0;@info>$i;$i++){
						$info{$info[$i]}=$i;
				}
				(defined $info{'SYMBOL'}) or die "ERROR:CSQ format not readed change script line 34\n";
		}elsif($_ =~ /^#/){$comment.=$_;
		}else{
				my @line=split(/\t/,);
				if( $chr ne $line[0]){
						if($chr ne ""){close OUTV;close OUTB;}
						elsif($chr eq "X"){last;}
						$chr=$line[0];
						push(@chr,$chr);
						print "doing $chr output berfore_liftover\n";
						open(OUTV,"|gzip -c >/Volumes/cancer/exac/before_liftover/$chr.vcf.gz");
						open(OUTB,">/Volumes/cancer/exac/before_liftover_bed/$chr.bed");
						print OUTV "$comment";
				}
				my ($ac,$an,$csq);
				if($line[7]=~/AC=([0-9,]+);.+;AN=(\d+);.+;CSQ=([^;]+);/){
						($ac,$an,$csq)=($1,$2,$3);
				}else{die "ERROR:cannot extract AC,AN,CAQ from $_\n";}
				my $focal=0;
				my @csq=split(/,/,$csq);
				for(my $i=0;@csq > $i;$i++){
						my @info=split(/\|/,$csq[$i]);
						if(defined $driver_genes{$info[$info{'SYMBOL'}]}){$focal++;}
				}
				if($focal>0){
#						if((length($line[3])>1)||(length($line[4]) >1)){next;}
						if($line[4]=~/,/){
								my $ref;
								my @alt=split(/,/,$line[4]);
								my @ac=split(/,/,$ac);
								my (@altout,@acout)=();
								for(my $i=0;@alt > $i;$i++){
										if(length($line[3])==1 and length($alt[$i])==1){
												push(@altout,$alt[$i]);
												push(@acout,$ac[$i]);
												$ref=$line[3];
										}elsif(length($line[3]) == length($alt[$i])){
												$ref=substr($line[3],0,1);
												my $altnuc=substr($alt[$i],0,1);
												if($ref eq $altnuc){die "ERROR:un expected line of ref&alt $_\n";}
												push(@altout,$altnuc);
												push(@acout,$ac[$i]);
										}
								}
								if(scalar(@altout) ==0){next;}
								my $alt=join(",",@altout);
								$ac=join(",",@acout);
								print OUTV join("\t",@line[0..2])."\t$ref\t$alt\t$line[5]\t$line[6]\t$ac\t$an\n";
								my $end=$line[1] +1;
								print OUTB "chr$line[0]\t$line[1]\t$end\t$line[0]:$line[1]\n";
						}elsif(length($line[3])==1 and length($line[4])==1 and $line[4] ne "-"){
								print OUTV join("\t",@line[0..6])."\t$ac\t$an\n";
								my $end=$line[1] +1;
								print OUTB "chr$line[0]\t$line[1]\t$end\t$line[0]:$line[1]\n";
						}
				}
		}
}
close OUTV;
close OUTB;
close VCF;

open (OUT,"|gzip -c >/Volumes/cancer/exac/exac_nontcga_liftovered.vcf.gz");
foreach my $chr (@chr){
		print "liftover $chr\n";
		my %remap = map{chomp;my @c=split("\t");$c[0]=~s/^chr//;($c[3],"$c[0]:$c[1]")}`~/vep/samtools/bin/liftOver /Volumes/cancer/exac/before_liftover_bed/$chr.bed ~/.vep/hg19ToHg38.over.chain /dev/stdout /dev/null 2>/dev/null`;
		print "printing out $chr \n";
		open (IN,"gunzip -c /Volumes/cancer/exac/before_liftover/$chr.vcf.gz|");
		while(<IN>){
				if($_=~/^#/){
						if($chr eq "1"){print OUT "$_";}
						next;
				}else{
						chomp;
						my @line=split(/\t/,);
						my ($chr,$pos,$refseq);
						if($remap{"$line[0]:$line[1]"} and $remap{"$line[0]:$line[1]"}=~ /^([^:]+):(\d+)$/){
								($chr,$pos)=($1,$2);
						}else{print  "WARNIG:there cant liftover $_\tso we dont use this variation\n";next;
						}
						my $fasta=`samtools faidx /Users/kazuki/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz $chr:$pos-$pos`;
						if($fasta =~ /\n([ATGC])\n$/){$refseq=$1;
						}else{die "ERROR: what happen?\n$_\n$fasta";}
						if($refseq eq $line[3]){
								if($line[4] =~ /,/){
										my @alt=split(/,/,$line[4]);
										my @ac =split(/,/,$line[7]);
										for(my $i=0;@alt>$i;$i++){
												print OUT "$chr\t$pos\t$line[2]\t$line[3]\t$alt[$i]\t$line[5]\t$line[6]\tAC=$ac[$i];AN=$line[8]\n";
										}
								}else{
										print OUT "$chr\t$pos\t".join("\t",@line[2..6])."\tAC=$line[7];AN=$line[8]\n";
								}
						}else{
								if($line[4] =~ /,/){
										my @alt=split(/,/,$line[4]);
										my @ac =split(/,/,$line[7]);
										for(my $i=0;@alt>$i;$i++){
												if($alt[$i] ne $refseq){
														print OUT "$chr\t$pos\t$line[2]\t$refseq\t$alt[$i]\t$line[5]\t$line[6]\tAC=$ac[$i];AN=$line[8]\n";
												}
										}
										my $refseq_ac = $line[8] - sum(@ac);
										print OUT "$chr\t$pos\t$line[2]\t$refseq\t$line[3]\t$line[5]\t$line[6]\tAC=$refseq_ac;AN=$line[8]\n";
								}else{
										if($line[4] ne $refseq){die "liftover correct?at $_\nliftover:".$remap{"$line[0]:$line[1]"}."\n";}
										my $refseq_ac = $line[8] - $line[7];
										print OUT "$chr\t$pos\t$line[2]\t$refseq\t$line[3]\t$line[5]\t$line[6]\tAC=$refseq_ac;AN=$line[8]\n";
								}
						}
				}
		}
		close IN;
}
close OUT;

exit;

										




