library(tidyr)
library(plyr)
library(dplyr)
library(pipeR)
library(stringr)
library(ggplot2)
library(readr)
library(readxl)
library(XML)
library(gtools)
library(purrr)
source("http://bioconductor.org/biocLite.R")
#biocLite()
library(GenomicRanges)
library(GenomicFeatures)

setwd('/working/maf/test/')
write_df_watal = function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  dplyr::select(gene,role)



mutect=read_tsv("mutect.maf.tsv",comment = "#") %>>%
  left_join(driver_genes,by=c("Hugo_Symbol"="gene")) %>>%
  filter(!is.na(role)) %>>%
  dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,
                Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Match_Norm_Seq_Allele1,
                Match_Norm_Seq_Allele2,HGVSc,HGVSp_Short,t_depth,t_ref_count,t_alt_count,n_depth,
                n_ref_count,n_alt_count,Consequence,CANONICAL,PolyPhen,FILTER)
varscan=read_tsv("./varscan.maf.tsv",comment = "#")%>>%
  left_join(driver_genes,by=c("Hugo_Symbol"="gene")) %>>%
  filter(!is.na(role)) %>>%
  dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,
                Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Match_Norm_Seq_Allele1,
                Match_Norm_Seq_Allele2,HGVSc,HGVSp_Short,t_depth,t_ref_count,t_alt_count,n_depth,
                n_ref_count,n_alt_count,Consequence,CANONICAL,PolyPhen,FILTER)
smsp=read_tsv("./somaticsniper.maf.tsv",comment = "#")%>>%
  left_join(driver_genes,by=c("Hugo_Symbol"="gene")) %>>%
  filter(!is.na(role)) %>>%
  dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,
                Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Match_Norm_Seq_Allele1,
                Match_Norm_Seq_Allele2,HGVSc,HGVSp_Short,t_depth,t_ref_count,t_alt_count,n_depth,
                n_ref_count,n_alt_count,Consequence,CANONICAL,PolyPhen,FILTER)
muse=read_tsv("./muse.maf.tsv",comment = "#") %>>%
  left_join(driver_genes,by=c("Hugo_Symbol"="gene")) %>>%
  filter(!is.na(role)) %>>%
  dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,
                Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Match_Norm_Seq_Allele1,
                Match_Norm_Seq_Allele2,HGVSc,HGVSp_Short,t_depth,t_ref_count,t_alt_count,n_depth,
                n_ref_count,n_alt_count,Consequence,CANONICAL,PolyPhen,FILTER)
norm_maf=read_tsv("/Volumes/areca42TB/tcga/maf_norm/brca/TCGA-5L-AAT1.maf",comment = "#") %>>%
  dplyr::select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,
                Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Match_Norm_Seq_Allele1,
                Match_Norm_Seq_Allele2,HGVSc,HGVSp_Short,t_depth,t_ref_count,t_alt_count,n_depth,
                n_ref_count,n_alt_count,Consequence,CANONICAL,PolyPhen,FILTER)

all = mutect %>>%
  filter(Variant_Type=="SNP") %>>%
  left_join(muse%>>%dplyr::select(Chromosome,Start_Position,FILTER)%>>%dplyr::rename(filter_muse=FILTER)) %>>%
  left_join(varscan%>>%dplyr::select(Chromosome,Start_Position,FILTER)%>>%dplyr::rename(filter_varscan=FILTER)) %>>%
  left_join(smsp%>>%dplyr::select(Chromosome,Start_Position,FILTER)%>>%dplyr::rename(filter_smsp=FILTER)) %>>%
  left_join(norm_maf%>>%dplyr::select(Chromosome,Start_Position)%>>%mutate(norm="have"))



#===========================================================================
#use when check 1kg maf and norm maf

classify_consequence = function(.data) {
  dplyr::mutate(.data, mutype= dplyr::recode(Consequence,
                                             downstream_gene_variant = 'flank',
                                             `3_prime_UTR_variant` = 'flank',
                                             upstream_gene_variant = 'flank',
                                             `5_prime_UTR_variant` = 'flank',
                                             frameshift_variant = 'truncating',
                                             frameshift_variant = 'truncating',
                                             inframe_deletion = 'inframe_indel',
                                             inframe_insertion = 'inframe_indel',
                                             intron_variant = 'silent',
                                             splice_region_variant = 'splice',
                                             coding_sequence_variant = 'missense',
                                             missense_variant = 'missense',
                                             stop_gained = 'truncating',
                                             stop_lost = 'truncating',
                                             stop_retained_variant = 'silent',
                                             synonymous_variant = 'silent',
                                             splice_acceptor_variant = 'splice',
                                             splice_donor_variant = 'splice',
                                             protein_altering_variant = 'missense',
                                             start_lost = 'truncating',
                                             `splice_region_variant,intron_variant` = 'splice',
                                             `stop_gained,frameshift_variant` = 'truncating',
                                             `splice_region_variant,synonymous_variant`='splice',
                                             `splice_region_variant,5_prime_UTR_variant`='splice',
                                             `missense_variant,splice_region_variant`='missense',
                                             `intron_variant,non_coding_transcript_variant`='silent',
                                             `non_coding_transcript_exon_variant,non_coding_transcript_variant`='silent'))
}

maf_root='/working/1000genomes/maf'
.infiles = list.files(maf_root, '.maf$', recursive=TRUE, full.names=TRUE) %>>% (?.)

maf_list = data_frame(path=.infiles, base=basename(path)) %>>%
  tidyr::separate(base, c('sample_allele', 'no'), '\\.', extra='drop') %>>%
  dplyr::select(-no) %>>%
  tidyr::separate(sample_allele,c("sample","allele"),'_',extra = 'drop') %>>%
  (?.)

.colnames = c('Hugo_Symbol','Chromosome','Start_Position','Consequence','PolyPhen')
.cols = .colnames %>>%
{setNames(c('c','c','d','c','c'), .)} %>>%
{do.call(readr::cols_only, as.list(.))}

strip_maf = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    filter(!mutype=='silent') %>>%
    filter(!mutype=='flank')
}
maf_1kg=maf_list %>>%
  mutate(maf=purrr::map(path,~strip_maf(.))) %>>%
  select(-path) %>>%
  unnest()


ascat_focal=read_tsv('/Volumes/areca42TB/tcga/CNA/brca/cel/error_of_annotate_ascat.txt') %>>%
  mutate(focal_ascat="no")
maf_focal=read_tsv('/Volumes/areca42TB/tcga/maf_norm/brca/list_of_perfect_maf.tsv') %>>%
  left_join(ascat_focal) %>>%
  mutate(focal=ifelse(is.na(focal_ascat),focal,"no")) %>>%dplyr::select(-focal_ascat)

maf_cancernorm = maf_focal %>>%
  mutate(filename=paste('/Volumes/areca42TB/tcga/maf_norm/brca/',patient_id,'.maf',sep="")) %>>%
  mutate(purrr::map(filename,~strip_maf(.))) %>>%
  unnest()
