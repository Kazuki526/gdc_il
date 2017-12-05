#!/usr/bin/perl
use warnings;
use strict;

my $single_tsv="/Volumes/areca42TB/GRCh38_singlefasta/GRCh38_single_proteincoding.tsv";
-e $single_tsv or die "ERROR::not existGRCh38_single_proteincoding.tsv\n";

my $gff="/Volumes/areca42TB/Homo_sapiens.GRCh38.84.gff3.gz";
-e $gff or die "ERROR::not exist GRCh38 gff3 file\n";


#pick single copy gene_symbol and ENST
open (IN,"$single_tsv");
my %singene=();
my %sinenst=();

