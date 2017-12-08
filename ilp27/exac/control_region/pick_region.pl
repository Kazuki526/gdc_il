#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);

########### perl pick_region.pl GENE_NAME.FILE OUT.DIR ######################

my $exac_path=$ENV{"HOME"}."/.vep/ExAC_nonTCGA.r1.sites.vep.vcf.gz";
(-e $exac_path) or die "ERROR:$exac_path is not exist!\n";
my $focal_genes_path=$ARGV[0];
(-e $focal_genes_path) or die "ERROR:$focal_genes_path is not exist!\n";
my $out_dir;
if(!defined $ARGV[1]){$out_dir=`pwd`;chomp $out_dir;
}else{$out_dir=$ARGV[1];}
(-e $out_dir) or die "ERROR::$out_dir is not exist!\n";

if(!-e "$out_dir/before_liftover"){
		print "there are not exist before_liftover dir so mkdir at $out_dir\n";
		mkdir "$out_dir/before_liftover";
}
if(!-e "$out_dir/before_liftover_indel"){
		mkdir "$out_dir/before_liftover_indel";
}
if(!-e "$out_dir/before_liftover_bed"){
		print "there are not exist before_liftover_bed dir so mkdir at $out_dir\n";
		mkdir "$out_dir/before_liftover_bed";
}
if(!-e "$out_dir/before_liftover_bed_indel"){
		mkdir "$out_dir/before_liftover_bed_indel";
}

#list up  focal_gene to %focal_genes
open(DG,"$focal_genes_path");
my %focal_genes=();
while(<DG>){
		if($_=~/^gene/){next;}
		chomp;
		my @line=split(/\t/,);
		if($line[0] eq "GDF2"){next;} ### <<<===この遺伝子はGRCh38で逆向き(逆転写?)にgenomeにあるのでややこしいからどかす!!
		$focal_genes{$line[0]}="ok";
}
close DG;

open(VCF,"gunzip -c $exac_path|");
my %csq_col_num=();
my @chr;
my ($chr,$comment)=("","");
my @info_ac=qw(AC AC_Adj AC_Het AC_Hom);
my @info_an=qw(AN AN_Adj);
foreach my $inf(qw(AC Het Hom)){
		foreach my $race(qw(AFR AMR EAS FIN NFE OTH SAS)){
				push(@info_ac,"$inf"."_$race");
		}
}
map{push(@info_an,"AN_".$_)}qw(AFR AMR EAS FIN NFE OTH SAS);

while(<VCF>){
		if($_ =~ /^##INFO=<ID=CSQ,.+Format:\s(.+)">$/){
				my @csq_col=split(/\|/,$1);
				for(my $i=0;@csq_col>$i;$i++){
						$csq_col_num{$csq_col[$i]}=$i;
				}
				(defined $csq_col_num{'SYMBOL'}) or die "ERROR:CSQ format not readed change script line 34\n";
		}elsif($_ =~ /^#/){$comment.=$_;
		}else{
				my @line=split(/\t/,);
				if( $chr ne $line[0]){
						if($chr ne ""){close OUTV;close OUTB;close OUTVI;close OUTBI;}
						elsif($chr eq "X"){last;}
						$chr=$line[0];
						push(@chr,$chr);
						print "doing $chr output berfore_liftover\n";
						open(OUTV,"|gzip -c >$out_dir/before_liftover/$chr.vcf.gz");
						open(OUTVI,"|gzip -c >$out_dir/before_liftover_indel/$chr.vcf.gz");
						open(OUTB,">$out_dir/before_liftover_bed/$chr.bed");
						open(OUTBI,">$out_dir/before_liftover_bed_indel/$chr.bed");
						print OUTV "$comment";print OUTVI "$comment";
				}		
				my %INFO=map{my@inf=split(/=/,);($inf[0],$inf[1])}split(/;/,$line[7]);
				my $focal=0;
				my @csq=split(/,/,$INFO{"CSQ"});
				foreach my $csq(@csq){
						my @csq_col=split(/\|/,$csq);
						if(defined $focal_genes{$csq_col[$csq_col_num{"SYMBOL"}]}){$focal++;}
				}
				if($focal>0){
						if($line[4]=~/,/){
								my $ref_snp="";
								my @alt=split(/,/,$line[4]);
								my @alt_snp=();
								my %ac_out=();
								for(my $i=0; @alt > $i; $i++){#普通のsnp
										if(length($line[3])==1 and length($alt[$i])==1){
												$ref_snp=$line[3];
												push(@alt_snp,$alt[$i]);
												foreach my $inf(@info_ac){
														if(defined $ac_out{$inf}){
																my @inf_ac=split(/,/,$INFO{$inf});
																$ac_out{$inf}.=",$inf_ac[$i]";
														}else{
																my @inf_ac=split(/,/,$INFO{$inf});
																$ac_out{$inf}=$inf_ac[$i];
														}
												}
										}elsif(length($line[3]) == length($alt[$i])){#ref=AGG,alt=TGGなどの時
												$ref_snp=substr($line[3],0,1);
												if($ref_snp eq substr($alt[$i],0,1)){die "ERROR:un expected line of ref&alt $_\n";}
														push(@alt_snp,substr($alt[$i],0,1));
														foreach my $inf(@info_ac){
																if(defined $ac_out{$inf}){
																		my @inf_ac=split(/,/,$INFO{$inf});
																		$ac_out{$inf}.=",$inf_ac[$i]";
																}else{
																		my @inf_ac=split(/,/,$INFO{$inf});
																		$ac_out{$inf}=$inf_ac[$i];
																}
														}
										}else{#indelの時
												my($posi,$ref_length,$alt_length,$ref,$alt)=
													($line[1],length($line[3]),length($alt[$i]),$line[3],$alt[$i]);
												while($ref_length>1 and $alt_length>1){
														($ref,$alt) = map{$_=substr($_,1)}($ref,$alt);
														--$ref_length;--$alt_length;++$posi;
												}
												print OUTVI join("\t",($line[0],$posi,".",$ref,$alt,$line[5],$line[6]))."\t";
												my @info_out=();
												map{my@inf_ac=split(/,/,$INFO{$_});push(@info_out,"$_=$inf_ac[$i]")}@info_ac;
												map{push(@info_out,"$_=$INFO{$_}")}@info_an;
												print OUTVI join(";",@info_out)."\n";
												my $end=$posi + $ref_length;
												print OUTBI "chr$line[0]\t$posi\t$end\t$line[0]:$posi-$end\n";
										}
								}
								if($ref_snp ne ""){
										print OUTV "$line[0]\t$line[1]\t.\t$ref_snp\t".join(",",@alt_snp)."\t$line[5]\t$line[6]\t";
										my @info_out=();
										map{push(@info_out,"$_=$ac_out{$_}")}@info_ac;
										map{push(@info_out,"$_=$INFO{$_}")}@info_an;
										print OUTV join(";",@info_out)."\n";
										my $end=$line[1] + 1;
										print OUTB "chr$line[0]\t$line[1]\t$end\t$line[0]:$line[1]-$end\n";
								}
						}else{
								my $info_out=join(";",map{"$_=$INFO{$_}"}(@info_ac,@info_an));
								my $end=$line[1] + length($line[3]);
								if(length($line[3])==1 and length($line[4])==1 and $line[4] ne "-"){
										print OUTV join("\t",@line[0..6])."\t$info_out\n";
										print OUTB "chr$line[0]\t$line[1]\t$end\t$line[0]:$line[1]-$end\n";
								}else{
										print OUTVI join("\t",@line[0..6])."\t$info_out\n";
										print OUTBI "chr$line[0]\t$line[1]\t$end\t$line[0]:$line[1]-$end\n";
								}
						}
				}
		}
}

close OUTV;close OUTB;close OUTVI;close OUTBI;
close VCF;


exit;
