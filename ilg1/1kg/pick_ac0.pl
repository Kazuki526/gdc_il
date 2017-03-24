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
print "reading bed file\n";
my %bed=();
open(BED,"$bed_file");
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0] =~ s/^chr//;
		for(my $i=$line[1];$line[2] >= $i;$i++){
				$bed{"$line[0]:$i"}="ok";
		}
}
close BED;

my @files=`ls $vcfpath|grep 'vcf.gz\$'`;chomp @files;
my @lines=();
foreach my $file (@files){
		push(@lines,`gunzip -c $vcfpath/$file|grep AC=0`);
		push(@lines,`gunzip -c $vcfpath/$file|grep 'AC=[0-9]\+,0'`);
		push(@lines,`gunzip -c $vcfpath/$file|grep 'AC=[0-9]\+,[0-9]\+,0;'`);
}
my $bed="/working/1000genomes/maf/extract/ac0.bed";
open(OUT,">$bed");
my %ac0=();
foreach my $line(@lines){
		chomp $line;
		my @line=split(/\t/,$line);
		my ($posi,$ref,$alt,$an);
		if($line[4]=~/,/){
				my @alt=split(/,/,$line[4]);
				my @ac;
				if($line[7] =~/AC=([0-9,]+);.+;AN=(\d+);/){
						@ac=split(/,/,$1);
						$an=$2;
				}else{die "ERROR:AC nodata at $line\n";}
				for(my $i=0;@alt>$i;$i++){
						if($ac[$i]==0){$alt=$alt[$i];}
				}
				if(length($line[3]) !=length($alt)){next;}
				elsif((length($line[3]) > 1)||(length($alt) > 1)){
						($posi,$ref,$alt)=&indel2snp($line[1],$line[3],$alt);
				}else{($posi,$ref)=($line[1],$line[3]);}
		}else{
				if($line[7]=~/AN=(\d+);/){
						$an=$1;
				}
				($posi,$ref,$alt)=($line[1],$line[3],$line[4]);
		}
		my $end=$posi+1;
		print OUT "chr$line[0]\t$posi\t$end\t$line[0]:$posi\n";
		$ac0{"$line[0]:$posi"}="$ref\t$alt\t0\t$an";
}
close OUT ;

my %remap= map{chomp;my @c=split(/\t/,);$c[0]=~s/^chr//;($c[3]=>"$c[0]:$c[1]")}`$liftover $bed $chainfile /dev/stdout /dev/null 2>/dev/null`;
open(OUT,">/working/1000genomes/maf/extract/ac0_variation.tsv");
print OUT "chr\tstart\tref\talt\tac_1kg\tan_1kg\n";
foreach my $region(sort(keys %ac0)){
		if(defined $bed{$remap{$region}}){
				my ($chr,$start)=split(/:/,$remap{$region});
				my ($ref,$alt,$ac,$an)=split(/\t/,$ac0{$region});
				my $faindex=`samtools faidx $hg38 $chr:$start-$start`;
				my ($fa_region,$ref_liftovered)=split(/\n/,$faindex,2);
				$ref_liftovered=~s/\s//;
				if($ref ne $ref_liftovered){
						if($alt ne $ref_liftovered){$ref=$ref_liftovered;
						}else{next;}
				}
				print OUT "$chr\t$start\t$ref\t$alt\t$ac\t$an\n";
		}
}
close OUT;


exit;



########################################################################
sub indel2snp {
		my ($start,$ref,$alt)=@_;
		while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
				($ref, $alt)=map{$_=substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $alt );
				$start++;
		}
		while( $ref and $alt and substr( $ref, -1, 1 ) eq substr( $alt, -1, 1 ) and $ref ne $alt ) {
				( $ref, $alt ) = map{$_ = substr( $_,0,-1  ); ( $_ ? $_ : "-" )} ( $ref, $alt );
		}
		return($start,$ref,$alt);
}

