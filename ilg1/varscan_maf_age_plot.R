library(tidyr)
library(plyr)
library(dplyr)
library(pipeR)
library(stringr)
library(ggplot2)
library(gridExtra)
library(readr)
library(readxl)
library(XML)
library(gtools)
library(purrr)

setwd('/Volumes/areca42TB/tcga/')
write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
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
topdriver_bed=read_tsv("/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed",col_names = c("chr","gene_start","gene_end","ids","score","strand"))%>>%
  tidyr::separate(ids,c("gene_symbol","ensg"),sep=";")
#brca_maf=read_tsv("/working/maf/7a07a833-4eab-44c9-bbf6-a64bd51c012e/TCGA.BRCA.muse.7a07a833-4eab-44c9-bbf6-a64bd51c012e.protected.maf.gz",comment = "#")
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
                                             protein_altering_variant = 'inframe_indel',
                                             start_lost = 'truncating',
                                             `splice_region_variant,intron_variant` = 'splice',
                                             `stop_gained,frameshift_variant` = 'truncating',
                                             `splice_region_variant,synonymous_variant`='splice',
                                             `splice_region_variant,5_prime_UTR_variant`='splice',
                                             `missense_variant,splice_region_variant`='missense',
                                             `intron_variant,non_coding_transcript_variant`='silent',
                                             `non_coding_transcript_exon_variant,non_coding_transcript_variant`='silent',
                                             `splice_region_variant,3_prime_UTR_variant`='flank',
                                             `stop_gained,splice_region_variant`='truncating',
                                             `stop_gained,protein_altering_variant`='truncating',
                                             `frameshift_variant,splice_region_variant`='truncating',
                                             `inframe_deletion,splice_region_variant`='splice',
                                             `splice_acceptor_variant,intron_variant`='splice',
                                             `splice_donor_variant,coding_sequence_variant,intron_variant`='splice',
                                             `stop_gained,inframe_deletion`='truncating',
                                             `stop_gained,inframe_insertion`='truncating'))
}
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  mutate(role=ifelse(role=="oncogene/TSG","TSG",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))

########################################################
################ germline mutation #####################
########################################################
extract_norm_maf=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2',"Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2",
                'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position",
                't_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c', 'c','c','c', 'c','c','c','c','c', 'd','d','d','d','d','d'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      classify_consequence()
  }
  #.bp is bodypart of it cancer
  maf=read_tsv(paste('/Volumes/areca42Tb2/gdc/varscan/',.bp,'/gender_age.tsv',sep="")) %>>%
    mutate(filename=paste('/Volumes/areca42Tb2/gdc/varscan/',.bp,"/maf/",patient_id,".maf",sep="")) %>>%
    mutate(purrr::map(filename,~strip_maf(.))) %>>%
    unnest() %>>%#(?.)%>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,t_allele1=Tumor_Seq_Allele1,t_allele2=Tumor_Seq_Allele2,
                  n_allele1=Match_Norm_Seq_Allele1,n_allele2=Match_Norm_Seq_Allele2) %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand)) %>>%
    mutate(soma_ro_germ=ifelse((t_allele1==n_allele1)&&(t_allele2==n_allele2),"germline","somatic"),
           t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
           n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
    mutate(LOH=ifelse((t_genotype=="homo")&&(n_genotype=="hetero"),"LOH","no"))
  
}

data.frame(body_part=c("breast","brain","lung","kidney","colorectal")) %>>%head(1)%>>%
  mutate(purrr::map(body_part,~extract_norm_maf(.))) %>>%
  unnest() %>>%
  write_df("~/Dropbox/install/tvz/body_part.maf.gz")
  
data.frame(cancer_type=c("hnsc","ov","prad","thca","ucec")) %>>%
  mutate(purrr::map(cancer_type,~extract_norm_maf(.))) %>>%
  unnest()%>>%
  write_df("~/Dropbox/install/tvz/cancer_type.maf.gz")
