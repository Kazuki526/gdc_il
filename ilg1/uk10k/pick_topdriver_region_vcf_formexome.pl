#!/usr/bin/perl
use warnings;
use strict;

my $pwd=`pwd`;chomp $pwd;
print "doing on $pwd\n";
if($pwd ne "/Volumes/areca42TB/ega/file"){die "ERROR:doing on wrong directory\n";}

my $bed_file="/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed";
(-e $bed_file) or die "ERROR:topdriver bed file:$bed_file is moved?\n";
print "reading bed file\n";
my %bed=();
my @chr=();
open(BED,"$bed_file");
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0] =~ s/^chr//;
		for(my $posi=$line[1];$line[2]>=$posi;$posi++){
				$bed{$line[0]}{$posi}="ok";
		}
		if(!grep(@chr,$line[0])){push(@chr,"$line[0]");}
}
close BED;

my @ls=`ls allvcf`;chomp @ls;
my $file_num=0;
my $header="";
my %all_vcf=();
foreach my $file(@ls){
		$file_num++;
		my ($chr)=("");
		open(VCF,"gunzip -c allvcf/$file|")or die "cant open $file\n";
		print "read and count $file\n";
		while(<VCF>){
				if($_=~/^##/){
						if($file_num==1){$header.="$_";}
						next;
				}
				chomp;
				my @line=split(/\t/,);
				if($line[0] eq "#CHROM"){
						if($file_num==1){
								$header.=join("\t",@line[0..8])."\n";
						}
						next;
				}
				if(!defined $bed{$line[0]}{$line[1]}){next;}
				my @alt=split(/,/,$line[4]);
				if($file_num==1){
						foreach my $alt (@alt){
								$all_vcf{$line[0]}{$line[1]}{$alt}{'all'}=join("\t",@line[0..6])."\t.\t$line[8]";
								$all_vcf{$line[0]}{$line[1]}{$alt}{'AN'}=0;
								$all_vcf{$line[0]}{$line[1]}{$alt}{'AC'}=0;
								$all_vcf{$line[0]}{$line[1]}{$alt}{'AA'}=0;#alt homo patient count
								$all_vcf{$line[0]}{$line[1]}{$alt}{'RA'}=0;#hetero patient coutn
						}
				}
				for(my $t=1;@alt>=$t;$t++){
						my ($ac,$an,$aa,$ra)=(0,0,0,0);
						for(my $i=9;@line>$i;$i++){
								my @gp=split(/:/,$line[$i]);
								if($gp[0] eq "./."){next;}
								$an+=2;
								if($gp[0] eq "$t/$t"){$ac+=2;$aa++;}
								elsif(($gp[0] =~ /$t\/[^$t]/)&&($gp[0] =~ /[^$t]\/$t/)){$ac++;$ra++;}
						}
						$all_vcf{$line[0]}{$line[1]}{$alt[$t-1]}{'AC'}+=$ac;
						$all_vcf{$line[0]}{$line[1]}{$alt[$t-1]}{'AA'}+=$aa;
						$all_vcf{$line[0]}{$line[1]}{$alt[$t-1]}{'RA'}+=$ra;
						$all_vcf{$line[0]}{$line[1]}{$alt[$t-1]}{'AN'}+=$an;
				}
		}
}
close VCF;

open(OUT1,"|gzip -c >all_sample_topdriver_region.vcf");
print OUT1 "$header";
open(OUT2,"|gzip -c >all_sample_topdriver_region_likevcf.tsv");
print OUT2 "chr\tstart\tref\talt\tuk_ac\tuk_an\tuk_althomo\tuk_hetero\n";
foreach my $chr(@chr){
		foreach my $posi (sort{$a <=> $b}keys(%{$all_vcf{$chr}})){
				foreach my $alt (keys %{$all_vcf{$chr}{$posi}}){
						my @line=split(/\t/,$all_vcf{$chr}{$posi}{$alt}{'all'});
						my ($ac,$aa,$ra,$an)=($all_vcf{$line[0]}{$line[1]}{$alt}{'AC'},
											  $all_vcf{$line[0]}{$line[1]}{$alt}{'AA'},
											  $all_vcf{$line[0]}{$line[1]}{$alt}{'RA'},
											  $all_vcf{$line[0]}{$line[1]}{$alt}{'AN'});
						print OUT1 "$chr\t$posi\t$line[2]\t$line[3]\t$alt\t$line[5]\t$line[6]\tAC=$ac;AN=$an;Althomo=$aa;hetero=$ra\t$line[8]\n";
						print OUT2 "$chr\t$posi\t$line[3]\t$alt\t$ac\t$an\t$aa\t$ra\n";
				}
		}
}
close OUT1;
close OUT2;


exit;


