#!/usr/bin/perl
use strict;
use warnings;

my $vcfpath="/working/1000genomes/topdrivers";
(-e $vcfpath) or die "ERROR:not exist $vcfpath\n";

my $liftover=$ENV{'HOME'}."/liftover/liftover";
(-e $liftover) or die "ERROR:liftover not exist at $liftover\n";
my $chainfile=$ENV{"HOME"}."/liftover/hg19ToHg38.over.chain";
(-e $chainfile) or die "ERROR:chain file not exist at $chainfile\n";
my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my $bed_file="/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed";
(-e $bed_file) or die "ERROR:topdriver bed file:$bed_file is moved?\n";


my @files=`ls $vcfpath|grep 'vcf.gz\$'`;chomp @files;

open(OUT,">/working/1000genomes/maf/extarct/ac0_variation.tsv");
print OUT "chr\tstart\tref\talt\tac_1kg\tan_1kg\n";
foreach my $file (@files){
		my @lines=`gunzip -c $vcfpath/$file|grep AC=0`;
		push(@lines,`gunzip -c $vcfpath/$file|grep 'AC=[0-9]\+,0'`);
		push(@lines,`gunzip -c $vcfpath/$file|grep 'AC=[0-9]\+,[0-9]\+,0;'`);
		foreach my $line(@lines){
				chomp $line;
				my @line=split(/\t/,$line);
				my ($posi,$ref,$alt,$an);
				if($line[4]=~/,/){
						my @alt=split(/,/$line[4]);
						my @ac;
						if($line[7] =~/AC=([0-9,]+);/){
								@ac=split(/,/,$1);
								for(my $i=0;@alt>$i;$i++){
										if($ac[$i]==0){$alt=$alt[$i];}
								}
								if(length($line[3]) !=length($alt)){next;}
								elsif((length($line[3]) > 1)||(length($alt) > 1)){
										($posi,$ref,$alt)=&indel2snp($line[1],$line[3],$alt);
								}else{($posi,$ref)=($line[1],$line[3]);}
						}
				}else{
				if($line[7]=~/AN=(\d+);/){
						$an=$1;
				}
				print OUT "$line[0]\t$line[1]\t$line[3]\t$alt\t0\t$an\n";
		}
}
close OUT ;
exit;



########################################################################
sub indel2snp {
		my ($start,$ref,$alt)=@_;
		while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
				($ref, $alt)=map{$_=substr(( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $alt );
				$start++;
		}
		while( $ref and $alt and substr( $ref, -1, 1 ) eq substr( $alt, -1, 1 ) and $ref ne $alt ) {
				( $ref, $alt ) = map{$_ = substr( $_,-1, 1 ); ( $_ ? $_ : "-" )} ( $ref, $alt );
		}
		return($start,$ref,$alt);
}

