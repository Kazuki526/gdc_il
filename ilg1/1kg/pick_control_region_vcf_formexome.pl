#!/usr/bin/perl
use warnings;
use strict;

my $pwd=`pwd`;chomp $pwd;
print "doing on $pwd\n";
if($pwd ne "/working/1000genomes"){die "ERROR:doing on wrong directory\n";}

my $liftover=$ENV{'HOME'}."/liftover/liftover";
(-e $liftover) or die "ERROR:liftover not exist at $liftover\n";
my $chainfile=$ENV{"HOME"}."/liftover/hg19ToHg38.over.chain";
(-e $chainfile) or die "ERROR:chain file not exist at $chainfile\n";
my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my $bed_file="/Volumes/areca42TB/GRCh38_singlefasta/control_genes_exon_with_splice_site.bed";
(-e $bed_file) or die "ERROR:topdriver bed file:$bed_file is moved?\n";
print "reading control genes bed file\n";
my %bed=();
my @chr=();
open(BED,"$bed_file");
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0] =~ s/^chr//;
		if(($line[3] =~ /^GPRIN2/) || ($line[3] =~ /^GDF2/)){next;}
		if(defined $bed{$line[0]}){
				$bed{$line[0]}.=":$line[1]-$line[2]";
		}else{
				$bed{$line[0]}="$line[1]-$line[2]";
		}
		if(!grep($_ eq $line[0],@chr)){push(@chr,"$line[0]");}
}
close BED;

=pod
#make to liftover bed
foreach my $chr (@chr){
		my $ls = "20130502/ALL.chr" . $chr . ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		if($chr eq "X"){$ls = "20130502/ALL.chr" . $chr . ".phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz";}
		open(VCF,"gunzip -c $ls|") or die "ERROR:cant gunzip -c $ls\n";
		open(OUT,">toLiftover/chr$chr.bed");
		print "writing toLiftover bed files\n";
		while(<VCF>){
				if($_=~/^#/){next;}
				chomp;
				my @line=split(/\t/,);
				my $end=$line[1]+1;
				#each position donot changed chromosome so not print out hg19:chr
				print OUT "chr$line[0]\t$line[1]\t$end\t$line[1]\n";
		}
		close OUT;
		close VCF
}
=cut

open(OUT1,"|gzip -c >control_genes/all_sample_control_region.vcf.gz");
open(OUT2,"|gzip -c >control_genes/all_sample_control_region_likevcf.tsv.gz");
print OUT2 "chr\tstart\tref\talt\t1kg_ac\t1kg_an\t1kg_althomo\t1kg_hetero\n";
my $file_num=0;
foreach my $chr(@chr){
		my %focal=();
		print "liftover chr$chr\n";
		my %remap= map{chomp;my @c=split(/\t/,);$c[0]=~s/^chr//;($c[3]=>$c[1])}`$liftover toLiftover/chr$chr.bed $chainfile /dev/stdout /dev/null 2>/dev/null`;
		foreach my $posi(keys %remap){
				foreach my $focal_region(split(/:/,$bed{$chr})){
						my ($start,$end)=split(/-/,$focal_region);
						if(($start <= $remap{$posi})&&($remap{$posi} <= $end)){
								$focal{$chr}{$posi}=$remap{$posi};
								last;
						}
				}
		}
		my $file = "20130502/ALL.chr" . $chr . ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
		if($chr eq "X"){$file = "20130502/ALL.chr" . $chr . ".phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz";}
		$file_num++;
		open(VCF,"gunzip -c $file|")or die "cant open $file\n";
		print "read and count $file\n";
		while(<VCF>){
				if($_=~/^##/){
						if($file_num==1){print OUT1 "$_";}
						next;
				}
				chomp;
				my @line=split(/\t/,);
				if($line[0] eq "#CHROM"){
						if($file_num==1){
								print OUT1 join("\t",@line[0..7])."\n";
						}
						next;
				}
				if(!defined $focal{$line[0]}{$line[1]}){next;}
				my $posi=$focal{$line[0]}{$line[1]};
				my @alt=split(/,/,$line[4]);
				my $fasta=`samtools faidx $hg38 $line[0]:$posi-$posi`;
				my ($region,$ref)=split(/\n/,$fasta);
				if(grep{length($_) >1}($line[3],@alt)){#indelの場合
						if($line[3] !~ /^$ref/){print "WARNING:chr$line[0]:$posi ref_colum is not same with hg38 reference\n";}
				}else{#indelでない場合
						if($line[3] ne $ref ){
								my $altn=0;
								for(my $t=1;@alt >= $t;$t++){
										if($alt[$t-1] eq $ref){
												$altn=$t;
												$alt[$t-1]=$line[3];
												$line[3]=$ref;
										}
								}
								if($altn==0){die "ERROR:chr$chr:$posi not have alt seq of correct ref\n";}
								for(my $i=9;@line > $i;$i++){
										eval "\$line[$i]=~ tr/0:$altn/$altn:0/";
								}
						}
				}
				for(my $t=1;@alt>=$t;$t++){
						if($alt[$t-1] =~/[^ATGC-]/){next;}
						my ($ac,$an,$aa,$ra)=(0,0,0,0);
						for(my $i=9;@line>$i;$i++){
								$an+=2;
								if($line[$i] eq "$t|$t"){$ac+=2;$aa++;}
								elsif(($line[$i] =~ /$t\|[^$t]/)||($line[$i] =~ /[^$t]\|$t/)){$ac++;$ra++;}
						}
						my ($posl,$refl,$altl)=($posi,$line[3],$alt[$t-1]);
						if(length($refl) != length($altl)){
								if(($refl !~ /^$altl/)&&($altl !~ /^$refl/)){die "ERROR::$line[0],$posl,$refl,$altl what indel ??";}
								while($refl and $altl and substr($refl,0,1) eq substr($altl,0,1) and $refl ne $altl){
										($refl,$altl) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $refl, $altl );
										$posl++;
								}
								while($refl and $altl and substr($refl,-1,1) eq substr($altl,-1,1) and $refl ne $altl){
#						if(grep{length($_)==2}($ref,$alt) and length($ref) != length($alt)){print "there are same tail indel $_\n";}
										($refl,$altl) = map{$_ = substr( $_, 0, -1 ); ( $_ ? $_ : "-" )} ( $refl, $altl );
								}
						}
						print OUT1 "$line[0]\t$posi\t$line[2]\t$line[3]\t$alt[$t-1]\t$line[5]\t$line[6]\tAC=$ac;AN=$an;Althomo=$aa;hetero=$ra\n";
						print OUT2 "$line[0]\t$posl\t$ref\t$altl\t$ac\t$an\t$aa\t$ra\n";
				}
		}
		close VCF;
}

close OUT1;
close OUT2;


exit;


