#!/usr/bin/perl
use warnings;
use strict;

my $pwd=`pwd`;chomp $pwd;
print "doing on $pwd\n";
if($pwd ne "/Volumes/areca42TB/ega/file"){die "ERROR:doing on wrong directory\n";}

my $liftover=$ENV{'HOME'}."/liftover/liftover";
(-e $liftover) or die "ERROR:liftover not exist at $liftover\n";
my $chainfile=$ENV{"HOME"}."/liftover/hg19ToHg38.over.chain";
(-e $chainfile) or die "ERROR:chain file not exist at $chainfile\n";
my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my $bed_file="/Volumes/areca42TB/GRCh38_singlefasta/control_genes_exon_with_splice_site.bed";
(-e $bed_file) or die "ERROR:topdriver bed file:$bed_file is moved?\n";
print "reading bed file\n";
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

my @ls=`ls allvcf`;chomp @ls;
#make to liftover bed
my $chr="";
=pod
open(VCF,"gunzip -c allvcf/$ls[0]|") or die "ERROR:cant gunzip-c allvcf/$ls[0]\n";
print "writing toLiftover bed files\n";
while(<VCF>){
		if($_=~/^#/){next;}
		chomp;
		my @line=split(/\t/,);
		if($chr ne "$line[0]"){
				if($chr ne ""){close OUT;}
				if(!grep($_ eq $line[0],@chr)){last;}
				$chr=$line[0];
				open(OUT,">toLiftover/chr$chr.bed");
		}
		my $end=$line[1]+1;
		#each position donot changed chromosome so not print out hg19:chr
		print OUT "chr$line[0]\t$line[1]\t$end\t$line[1]\n";
}
=cut

close VCF;
my %focal=();
foreach $chr(@chr){
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
}


my $file_num=0;
my $header="";
my %all_vcf=();
foreach my $file(@ls){
		$file_num++;
		open(VCF,"gunzip -c allvcf/$file|")or die "cant open $file\n";
		print "read and count $file\n";
		my $fline=0;
		while(<VCF>){
				if($_=~/^##/){
						if($file_num==1){$header.="$_";}
						next;
				}
				chomp;
				my @line=split(/\t/,);
				if($line[0] eq "#CHROM"){
						if($file_num==1){
								$header.=join("\t",@line[0..7])."\n";
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
								if($altn==0){die "ERROR:chr$line[0]:$posi not have alt seq of correct ref\n";}
								for(my $i=9;@line > $i;$i++){
										eval "\$line[$i]=~ tr/0:$altn/$altn:0/";
								}
						}
				}
				if($file_num==1){
						foreach my $alt (@alt){
								$all_vcf{$line[0]}{$posi}{$alt}{'all'}=join("\t",@line[0..6])."\t.";
								$all_vcf{$line[0]}{$posi}{$alt}{'AN'}=0;
								$all_vcf{$line[0]}{$posi}{$alt}{'AC'}=0;
								$all_vcf{$line[0]}{$posi}{$alt}{'AA'}=0;#alt homo patient count
								$all_vcf{$line[0]}{$posi}{$alt}{'RA'}=0;#hetero patient coutn
						}
				}
				for(my $t=1;@alt>=$t;$t++){
						my ($ac,$an,$aa,$ra)=(0,0,0,0);
						for(my $i=9;@line>$i;$i++){
								my @gp=split(/:/,$line[$i]);
								if($gp[0] eq "./."){next;}
								$an+=2;
								if($gp[0] eq "$t/$t"){$ac+=2;$aa++;}
								elsif(($gp[0] =~ /$t\/[^$t]/)||($gp[0] =~ /[^$t]\/$t/)){$ac++;$ra++;}
						}
						$all_vcf{$line[0]}{$posi}{$alt[$t-1]}{'AC'}+=$ac;
						$all_vcf{$line[0]}{$posi}{$alt[$t-1]}{'AA'}+=$aa;
						$all_vcf{$line[0]}{$posi}{$alt[$t-1]}{'RA'}+=$ra;
						$all_vcf{$line[0]}{$posi}{$alt[$t-1]}{'AN'}+=$an;
				}
		}
}
close VCF;

open(OUT1,"|gzip -c >all_sample_control_region.vcf.gz");
print OUT1 "$header";
open(OUT2,"|gzip -c >all_sample_control_region_likevcf.tsv.gz");
print OUT2 "chr\tstart\tref\talt\tuk_ac\tuk_an\tuk_althomo\tuk_hetero\n";
foreach my $chr(@chr){
		foreach my $posi (sort{$a <=> $b}keys(%{$all_vcf{$chr}})){
				foreach my $alt (keys %{$all_vcf{$chr}{$posi}}){
						my @line=split(/\t/,$all_vcf{$chr}{$posi}{$alt}{'all'});
						my ($refv,$altv) = ($line[3],$alt);
						if((length($refv) != length($altv)) && (($refv !~/^$altv/)||($altv!~/^$refv/))){
								while($refv and $altv and substr($refv,-1,1) eq substr($altv,-1,1) and (($refv =~/^$altv/) || ($altv =~/^$refv/))){
										($refv,$altv) = map{$_ = substr( $_, 0, -1 ); ( $_ ? $_ : "-" )} ( $refv, $altv );
								}
						}
						my ($ac,$aa,$ra,$an)=($all_vcf{$line[0]}{$posi}{$alt}{'AC'},
											  $all_vcf{$line[0]}{$posi}{$alt}{'AA'},
											  $all_vcf{$line[0]}{$posi}{$alt}{'RA'},
											  $all_vcf{$line[0]}{$posi}{$alt}{'AN'});
						if($an==0){next;}
						my ($refl,$altl,$posil)=($line[3],$alt,$posi);
						if(length($refl) != length($altl)){
								while($refl and $altl and substr($refl,0,1) eq substr($altl,0,1) and $refl ne $altl){
										($refl,$altl) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $refl, $altl );
										$posil++;
								}
								while($refl and $altl and substr($refl,-1,1) eq substr($altl,-1,1) and $refl ne $altl){
#						if(grep{length($_)==2}($ref,$alt) and length($ref) != length($alt)){print "there are same tail indel $_\n";}
										($refl,$altl) = map{$_ = substr( $_, 0, -1 ); ( $_ ? $_ : "-" )} ( $refl, $altl );
								}
						}
						print OUT1 "$chr\t$posi\t$line[2]\t$refv\t$altv\t$line[5]\t$line[6]\tAC=$ac;AN=$an;Althomo=$aa;hetero=$ra\n";
						print OUT2 "$chr\t$posil\t$refl\t$altl\t$ac\t$an\t$aa\t$ra\n";
				}
		}
}
close OUT1;
close OUT2;


exit;


