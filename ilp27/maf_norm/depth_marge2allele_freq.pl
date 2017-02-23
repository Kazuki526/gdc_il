#!/usr/bin/perl
use warnings;
use strict;

#perl depth_by_bam.pl body_part 

my $pwd=`pwd`;chomp $pwd;
if($pwd ne "/Volumes/cancer/kaz_gdc"){die "!!ERROR!! doing on wrong dir\n";}

my $bp=$ARGV[0]; #body part

my @ls=`ls $bp/depth/`;
my %af=();
foreach my $file(@ls){
		chomp $file;
		if($file eq "all_coverage.tsv"){next;}
		print "input $file\n";
		open(IN,"$bp/depth/$file");
		<IN>;  #delete colum line
		while(<IN>){
				chomp;
				my @line=split(/\t/,);
				$line[0] =~ s/^chr//;
				if(!defined $af{$line[0]}{$line[1]}){$af{$line[0]}{$line[1]}=0;}
				for(my $i=2;@line>$i;$i++){
						if ($line[$i] > 4){$af{$line[0]}{$line[1]}++;}
				}
		}
		close IN;
}

open(OUT,">$bp/depth/all_coverage.tsv");
print OUT "Chromosome\tPosition\tcoverage\n";
foreach my $chr (sort{$a<=>$b} keys %af){
		foreach my $posi(sort{$a<=>$b} keys %{$af{$chr}}){
				print OUT "$chr\t$posi\t$af{$chr}{$posi}\n";
		}
}
close OUT;

exit;

