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
  mutate(role=ifelse(role=="oncogene/TSG","oncogene",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))

if(0){  ##################  from here brefore script
.colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
              'Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count','Consequence','PolyPhen')
.cols = .colnames %>>%
{setNames(c('c','c','d','d','c','c','c','d','d','d','c','c'), .)} %>>%
{do.call(readr::cols_only, as.list(.))}
strip_maf = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    filter(!mutype=='silent') %>>%
    filter(!mutype=='flank')
}

norm_maf=maf_focal %>>%
  filter(focal=="ok") %>>% #head(100) %>>%
  mutate(filename=paste('maf_norm/breast/',patient_id,'.maf',sep="")) %>>%
  mutate(purrr::map(filename,~strip_maf(.))) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                ref=Reference_Allele,t_alle1=Tumor_Seq_Allele1,t_alle2=Tumor_Seq_Allele2) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  dplyr::select(patient_id,gene_symbol,chr,start,end,ref,t_alle1,t_alle2,mutype,PolyPhen,
                t_depth,t_ref_count,t_alt_count)

norm_maf %>>%
  #dplyr::select(t_depth,t_ref_count,t_alt_count) %>>%
  mutate(t_alt2_count=t_depth - t_ref_count - t_alt_count) %>>%
  View

norm_tally = norm_maf %>>%
  tidyr::gather(allele,alt,t_alle1,t_alle2) %>>%
  filter(ref != alt) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen) %>>%
  summarise(n_norm=n())
} ##########################################to here

########################################################
################ germline mutation #####################
########################################################
extract_norm_maf=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2','Consequence','PolyPhen')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      classify_consequence() %>>%
      filter(!mutype=='silent') %>>%
      filter(!mutype=='flank')
  }
#.bp is bodypart of it cancer
  maf=read_tsv(paste('maf_norm/',.bp,'/depth/gender_file.tsv',sep="")) %>>%
    mutate(filename=paste('maf_norm/',.bp,"/",patient_id,".maf",sep="")) %>>%
    mutate(purrr::map(filename,~strip_maf(.))) %>>%
    unnest() %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,allele1=Tumor_Seq_Allele1,allele2=Tumor_Seq_Allele2) %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand)) %>>%
    tidyr::gather(allele,alt,allele1,allele2) %>>%
    filter(!ref==alt)%>>%
    filter(!(gender=="male"&chr=="X"&allele=="allele2"))%>>%
    dplyr::select(patient_id,gene_symbol,allele,chr,start,end,ref,alt,mutype,PolyPhen)
  
}
#maf_norm_breast=extract_norm_maf("breast")
#write_df(maf_norm_breast,"maf_norm/united_maf/breast_nonsilent.maf")
maf_norm_breast=read_tsv("maf_norm/united_maf/breast_nonsilent.maf")

#maf_norm_lung=extract_norm_maf("lung")
#write_df(maf_norm_lung,"maf_norm/united_maf/lung_nonsilent.maf")
maf_norm_lung=read_tsv("maf_norm/united_maf/lung_nonsilent.maf")

#maf_norm_brain=extract_norm_maf("brain")
#write_df(maf_norm_brain,"maf_norm/united_maf/brain_nonsilent.maf")
maf_norm_brain=read_tsv("maf_norm/united_maf/brain_nonsilent.maf")

#maf_norm_kidney=extract_norm_maf("kidney")
#write_df(maf_norm_kidney,"maf_norm/united_maf/kidney_nonsilent.maf")
maf_norm_kidney=read_tsv("maf_norm/united_maf/kidney_nonsilent.maf")

#maf_norm_colorectal=extract_norm_maf("colorectal")
#write_df(maf_norm_colorectal,"maf_norm/united_maf/colorectal_nonsilent.maf")
maf_norm_colorectal=read_tsv("maf_norm/united_maf/colorectal_nonsilent.maf")

#########silent###############
extract_norm_maf_silent=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2','Consequence','PolyPhen')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      filter(Consequence=="synonymous_variant")
  }
  #.bp is bodypart of it cancer
  maf=read_tsv(paste('maf_norm/',.bp,'/depth/gender_file.tsv',sep="")) %>>%
    mutate(filename=paste('maf_norm/',.bp,"/",patient_id,".maf",sep="")) %>>%
    mutate(purrr::map(filename,~strip_maf(.))) %>>%
    unnest() %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,allele1=Tumor_Seq_Allele1,allele2=Tumor_Seq_Allele2) %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand)) %>>%
    tidyr::gather(allele,alt,allele1,allele2) %>>%
    filter(!ref==alt)%>>%
    filter(!(gender=="male"&chr=="X"&allele=="allele2"))%>>%
    dplyr::select(patient_id,gene_symbol,allele,chr,start,end,ref,alt,PolyPhen)
  
}
#maf_norm_breast_silent=extract_norm_maf_silent("breast")
#write_df(maf_norm_breast_silent,"maf_norm/united_maf/breast_silent.maf")
maf_norm_breast_silent=read_tsv("maf_norm/united_maf/breast_silent.maf")

#maf_norm_lung_silent=extract_norm_maf_silent("lung")
#write_df(maf_norm_lung_silent,"maf_norm/united_maf/lung_silent.maf")
maf_norm_lung_silent=read_tsv("maf_norm/united_maf/lung_silent.maf")

#maf_norm_brain_silent=extract_norm_maf_silent("brain")
#write_df(maf_norm_brain_silent,"maf_norm/united_maf/brain_silent.maf")
maf_norm_brain_silent=read_tsv("maf_norm/united_maf/brain_silent.maf")

#maf_norm_kidney_silent=extract_norm_maf_silent("kidney")
#write_df(maf_norm_kidney_silent,"maf_norm/united_maf/kidney_silent.maf")
maf_norm_kidney_silent=read_tsv("maf_norm/united_maf/kidney_silent.maf")

#maf_norm_colorectal_silent=extract_norm_maf_silent("colorectal")
#write_df(maf_norm_colorectal_silent,"maf_norm/united_maf/colorectal_silent.maf")
maf_norm_colorectal_silent=read_tsv("maf_norm/united_maf/colorectal_silent.maf")
###########################################
######### somatic mutaion #################
###########################################

#if .bp=lung it is lusc luad is correct
#       brain it is lgg gbm is correct
#       kidney it is kirc kirp kich
muse_mutect_every_path=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2','Consequence','PolyPhen','FILTER','Tumor_Sample_Barcode',
                't_depth','t_ref_count','t_alt_count','CDS_position','HGVSc','Protein_position','HGVSp_Short',
                'CANONICAL','Transcript_ID')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c','c','c','c','c','c','d','d','d','c','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      classify_consequence() %>>%
      filter(FILTER=="PASS") %>>%
      dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                    ref=Reference_Allele,allele1=Tumor_Seq_Allele1,allele2=Tumor_Seq_Allele2) %>>%
      left_join(topdriver_bed %>>% dplyr::select(gene_symbol,strand)) %>>%
      filter(!is.na(strand)) %>>%
      mutate(patient_id=str_extract(Tumor_Sample_Barcode,"TCGA-[^-]*-[^-]*"))%>>%
      filter(!mutype=='silent') %>>%
      filter(!mutype=='flank')
  }
  .BP=toupper(.bp)
  .maf=read_tsv("/working/maf/GDC_controlled_MAFs_manifest.txt") %>>%
    filter(str_detect(filename,.BP)) %>>%
    filter(str_detect(filename,"mutect")|str_detect(filename,"muse")) %>>%
    mutate(path=paste("/working/maf",id,filename,sep="/")) %>>%
    mutate(maf=purrr::map(path,~strip_maf(.)))
  out_maf=full_join(.maf$maf[[1]],.maf$maf[[2]])
}

#lusc_maf=muse_mutect_every_path("lusc")%>>%mutate(cancer_type="LUSC")
#write_df(lusc_maf,"kaz_maf/extracted_maf/lusc_topdriver105genes.maf")
lusc_maf=read_tsv("kaz_maf/extracted_maf/lusc_topdriver105genes.maf")
#luad_maf=muse_mutect_every_path("luad")%>>%mutate(cancer_type="LUAD")
#write_df(luad_maf,"kaz_maf/extracted_maf/luad_topdriver105genes.maf")
luad_maf=read_tsv("kaz_maf/extracted_maf/luad_topdriver105genes.maf")
lung_maf=full_join(lusc_maf,luad_maf) %>>%
  left_join(read_tsv("maf_norm/lung/depth/gender_file.tsv")) %>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)

#breast_maf=muse_mutect_every_path("brca") %>>%mutate(cancer_type="BRCA")
#write_df(breast_maf,"kaz_maf/extracted_maf/brca_topdriver105genes.maf")
breast_maf=read_tsv("kaz_maf/extracted_maf/brca_topdriver105genes.maf")

#lgg_maf=muse_mutect_every_path("lgg")%>>%mutate(cancer_type="LGG")
#write_df(lgg_maf,"kaz_maf/extracted_maf/lgg_topdriver105genes.maf")
lgg_maf=read_tsv("kaz_maf/extracted_maf/lgg_topdriver105genes.maf")
#gbm_maf=muse_mutect_every_path("gbm")%>>%mutate(cancer_type="GBM")
#write_df(gbm_maf,"kaz_maf/extracted_maf/gbm_topdriver105genes.maf")
gbm_maf=read_tsv("kaz_maf/extracted_maf/gbm_topdriver105genes.maf")
brain_maf=full_join(lgg_maf,gbm_maf) %>>%
  left_join(read_tsv("maf_norm/brain/depth/gender_file.tsv")) %>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)

#muse_mutect_every_path("kirc")%>>%mutate(cancer_type="KIRC") %>>%
#  write_df("kaz_maf/extracted_maf/kirc_topdriver105genes.maf")
kirc_maf=read_tsv("kaz_maf/extracted_maf/kirc_topdriver105genes.maf")
#muse_mutect_every_path("kirp")%>>%mutate(cancer_type="KIRP") %>>%
#  write_df("kaz_maf/extracted_maf/kirp_topdriver105genes.maf")
kirp_maf=read_tsv("kaz_maf/extracted_maf/kirp_topdriver105genes.maf")
muse_mutect_every_path("kich")%>>%mutate(cancer_type="KICH") %>>%
  write_df("kaz_maf/extracted_maf/kich_topdriver105genes.maf")
kich_maf=read_tsv("kaz_maf/extracted_maf/kich_topdriver105genes.maf")
kidney_maf=full_join(kirc_maf,kirp_maf)%>>%full_join(kich_maf) %>>%
  left_join(read_tsv(("maf_norm/kidney/depth/gender_file.tsv")))%>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)

muse_mutect_every_path("coad")%>>%mutate(cancer_type="COAD") %>>%
  write_df("kaz_maf/extracted_maf/coad_topdriver105genes.maf")
coad_maf=read_tsv("kaz_maf/extracted_maf/coad_topdriver105genes.maf")
muse_mutect_every_path("read")%>>%mutate(cancer_type="READ") %>>%
  write_df("kaz_maf/extracted_maf/read_topdriver105genes.maf")
read_maf=read_tsv("kaz_maf/extracted_maf/read_topdriver105genes.maf")
colorectal_maf=full_join(coad_maf,read_maf)%>>%
  left_join(read_tsv(("maf_norm/colorectal/depth/gender_file.tsv")))%>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)

muse_mutect_every_path('ov')%>>%mutate(cancer_type='OV') %>>%
  write_df("kaz_maf/extracted_maf/ov_topdriver105genes.maf")
muse_mutect_every_path('hnsc')%>>%mutate(cancer_type='HNSC') %>>%
  write_df("kaz_maf/extracted_maf/hnsc_topdriver105genes.maf")
muse_mutect_every_path('prad')%>>%mutate(cancer_type='PRAD') %>>%
  write_df("kaz_maf/extracted_maf/prad_topdriver105genes.maf")
muse_mutect_every_path('thca')%>>%mutate(cancer_type='THCA') %>>%
  write_df("kaz_maf/extracted_maf/thca_topdriver105genes.maf")
muse_mutect_every_path('ucec')%>>%mutate(cancer_type='UCEC') %>>%
  write_df("kaz_maf/extracted_maf/ucec_topdriver105genes.maf")

###########################################
################  1kg  ####################
###########################################
if(0){#needed to do only first time
maf_root='/working/1000genomes/maf'
.infiles = list.files(maf_root, '.maf$', recursive=TRUE, full.names=TRUE) %>>% (?.)

maf_list = data_frame(path=.infiles, base=basename(path)) %>>%
  tidyr::separate(base, c('sample_allele', 'no'), '\\.', extra='drop') %>>%
  dplyr::select(-no) %>>%
  tidyr::separate(sample_allele,c("sample","allele"),'_',extra = 'drop') %>>%
  (?.)

.colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position',
              'Reference_Allele','Tumor_Seq_Allele2','Consequence','PolyPhen')
.cols = .colnames %>>%
{setNames(c('c','c','d','d','c','c','c','c'), .)} %>>%
{do.call(readr::cols_only, as.list(.))}

strip_maf = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    filter(!mutype=='silent') %>>%
    filter(!mutype=='flank')
}

have_y_sample=read_tsv("/working/1000genomes/topdrivers/have_chry_samples.tsv") %>>%
  mutate(allele="allele2",chr="X")

maf_1kg_nonsilent=maf_list %>>%
  mutate(maf=purrr::map(path,~strip_maf(.))) %>>%
  select(-path) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,
                end=End_Position,ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  dplyr::select(sample,allele,gene_symbol,chr,start,end,ref,alt,mutype,PolyPhen) %>>%
  left_join(have_y_sample) %>>%
  filter(is.na(Y)) %>>%
  dplyr::select(-Y)
write_df(maf_1kg_nonsilent,"/working/1000genomes/maf/extract/nonsilent.maf")
}
maf_1kg_nonsilent=read_tsv("/working/1000genomes/maf/extract/nonsilent.maf")

vcf_1kg_ac0=read_tsv("/working/1000genomes/maf/extract/ac0_variation.tsv")

#####silent#####
strip_maf_silent = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    filter(Consequence=="synonymous_variant")
}
#maf_1kg_silent=maf_list %>>%
  mutate(maf=purrr::map(path,~strip_maf_silent(.))) %>>%
  select(-path) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,
                end=End_Position,ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
  dplyr::select(sample,allele,gene_symbol,chr,start,end,ref,alt,PolyPhen) %>>%
  left_join(have_y_sample) %>>%
  filter(is.na(Y)) %>>%
  dplyr::select(-Y)

#write_df(maf_1kg_silent,"/working/1000genomes/maf/extract/silent.maf")
maf_1kg_silent=read_tsv("/working/1000genomes/maf/extract/silent.maf") 

##############################################
###########UK 10K ############################
##############################################
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_likevcf.tsv.gz",col_types = "cdccdddd")


##############################################
############# EXAC nonTCGA ###################
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",col_types = "cdccdd")

###################################################################################################################
################## see all bam coverage###################

lung_depth_byexon=read_tsv("maf_norm/lung/depth/exon_mean_depth_tidy.tsv.gz")
lung_depth_byexon_ =lung_depth_byexon %>>%filter(patient_id=="TCGA-05-4244") %>>%mutate(exon_length=end - start)

##############################################################################################################
################################################  testing   ##################################################
##############################################################################################################
###fisher exact test by site
sample_num_1kg=5008
sample_num_X_1kg=sample_num_1kg - 1233

lung_coverage =read_tsv("/Volumes/areca42TB/tcga/maf_norm/lung/depth/coverage_all_data_exist_patient.tsv",
                       col_names = c("chr","start","an_lung"),col_types = "cdd")
brain_coverage=read_tsv("/Volumes/areca42TB/tcga/maf_norm/brain/depth/coverage_all_data_exist_patient.tsv",
                        col_names = c("chr","start","an_brain"),col_types = "cdd")
tally_lung=maf_norm_lung %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen) %>>%
  summarise(ac_lung=n()) 

tally_1kg = maf_1kg %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen) %>>%
  summarise(ac_1kg=n()) %>>%
  full_join(vcf_1kg_ac0)

test_fisher =function(.row){
  .test=fisher.test(matrix(c(.row$ac_1kg,.row$ac_lung,max(1,.row$an_1kg-.row$ac_1kg),max(1,.row$an_lung-.row$ac_lung)),2))
  .test$p.value
}


variatn_lung=full_join(tally_1kg,tally_lung) %>>%#head(10)%>>%
  mutate(ac_1kg=as.double(ifelse(is.na(ac_1kg),0,ac_1kg)),ac_lung=as.double(ifelse(is.na(ac_lung),0,ac_lung))) %>>%
  left_join(lung_coverage) %>>%
  filter(!is.na(an_lung)) %>>%
  #filter(ac_lung>2,ac_1kg>2) %>>%
  mutate(an_1kg=ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg)) %>>%
  by_row(.to="p_value",~test_fisher(.)) %>>% unnest() %>>%
  mutate(compar=(ac_lung / an_lung)/(ac_1kg / an_1kg))

variatn_lung %>>%filter(ac_lung > 2,ac_1kg > 2) %>>%
  ggplot()+
  geom_histogram(aes(x=p_value,fill=mutype),bins=30)+
  scale_x_log10()

#######silent#######
tally_lung_silent=maf_norm_lung_silent %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,PolyPhen) %>>%
  summarise(ac_lung=n()) 

tally_1kg_silent = maf_1kg_silent %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,PolyPhen) %>>%
  summarise(ac_1kg=n())



all_variant_lung_silent=full_join(tally_1kg_silent,tally_lung_silent) %>>%#head(10)%>>%
  mutate(ac_1kg=as.double(ifelse(is.na(ac_1kg),0,ac_1kg)),ac_lung=as.double(ifelse(is.na(ac_lung),0,ac_lung))) %>>%
  left_join(lung_coverage) %>>%
  filter(!is.na(an_lung)) %>>%
  #filter(ac_lung>2,ac_1kg>2) %>>%
  mutate(an_1kg=ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg)) %>>%
  by_row(.to="p_value",~test_fisher(.)) %>>% unnest() %>>%
  mutate(compar=(ac_lung / an_lung)/(ac_1kg / an_1kg))

all_variant_lung_silent %>>%filter(ac_lung > 2,ac_1kg > 2) %>>%
  filter(compar<10,compar>0.1) %>>%
  ggplot()+
  geom_histogram(aes(x=p_value),binwidth = 0.01)


######################## 有意差のあったsiteを持っている人はmutationが少ない？#################################
#上の all_variant_~~はそのまま使用
#lung_maf,maf_norm_lung

#CNA でLOHのものをtruncatingとできる。
lung_asact=read_tsv("CNA/lung/cel/annotate_ascat.tsv.gz") %>>%
  group_by(patient_id,gene_symbol) %>>%
  summarise(nminor=min(nminor),nmajor=min(nmajor)) %>>%
  filter(nminor==0)
by_site=variatn_lung %>>%#head(10)%>>%
  #filter(p_valu<0.0001) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  by_row(function(.row){
    .to_fisher=read_tsv("maf_norm/lung/depth/gender_file.tsv") %>>%
      left_join(maf_norm_lung %>>%
                  filter(chr==.row$chr,start==.row$start,end==.row$end,alt==.row$alt)%>>%
                  dplyr::select(patient_id) %>>% mutate(germline="have")) %>>%
      left_join(full_join(lung_maf %>>%filter(gene_symbol==.row$gene_symbol)%>>%
                  #filter(mutype=="missense")%>>%
                  dplyr::select(patient_id)%>>%
                  dplyr::distinct(patient_id) %>>% mutate(somatic="have"),
                lung_asact %>>%filter(gene_symbol==.row$gene_symbol,nmajor == 1)%>>%
                  dplyr::select(patient_id)%>>%mutate(somatic="have"))
                  )%>>%
      mutate(somatic=ifelse(is.na(somatic),"no","have"),germline=ifelse(is.na(germline),"no","have")) %>>%
      group_by(germline,somatic)%>>%tally() %>>%
      dplyr::arrange(germline,somatic)#これでghsh,ghsn,gnsh,gnsnの順になる(g=germ,s=somatic,h=have,n=no)
    if(length(.to_fisher$n) < 4){
      Inf
    }else{
      .fisher=fisher.test(matrix(c(.to_fisher$n),2))
      .fisher$p.value
    }
  },.to="by_site_pvalue") %>>%
  unnest()

by_site %>>%filter(by_site_pvalue!=Inf)%>>%
  filter(p_value < 0.1) %>>%
  ggplot(aes(x=p_valu,y=by_site_pvalue))+geom_point()+scale_x_log10()#+scale_y_log10()

by_site %>>%
  filter(by_site_pvalue != Inf) %>>%
  #filter(p_value < 0.1) %>>%
  {cor.test(.$p_value,.$by_site_pvalue)}

#######　上の検定をmutationのあるgeneだけでなく全遺伝子で###########
by_sitex105genes=variatn_lung %>>%
  filter(ac_lung>2,ac_1kg>2) %>>%
  #filter(p_valu<0.0001) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%#head(10)%>>%
  by_row(function(.row){
    driver_genes %>>%
    dplyr::select(gene,role) %>>%
    by_row(function(.gene){
      .to_fisher=read_tsv("maf_norm/lung/depth/gender_file.tsv") %>>%
        left_join(maf_norm_lung %>>%
                    filter(chr==.row$chr,start==.row$start,end==.row$end,alt==.row$alt)%>>%
                    dplyr::select(patient_id) %>>% mutate(germline="have")) %>>%
        left_join(full_join(lung_maf %>>%filter(gene_symbol==.gene$gene)%>>%
                              #filter(mutype=="missense")%>>%
                              dplyr::select(patient_id)%>>%
                              dplyr::distinct(patient_id) %>>% mutate(somatic="have"),
                            lung_asact %>>%filter(gene_symbol==.gene$gene,nmajor == 1)%>>%
                              dplyr::select(patient_id)%>>%mutate(somatic="have"))
        )%>>%
        mutate(somatic=ifelse(is.na(somatic),"no","have"),germline=ifelse(is.na(germline),"no","have")) %>>%
        group_by(germline,somatic)%>>%tally() %>>%
        dplyr::arrange(germline,somatic)#これでghsh,ghsn,gnsh,gnsnの順になる(g=germ,s=somatic,h=have,n=no)
      if(length(.to_fisher$n) < 4){
        Inf
      }else{
        .fisher=fisher.test(matrix(c(.to_fisher$n),2))
        .fisher$p.value
      }
    },.to="genes_pvalue")%>>%unnest()%>>%filter(genes_pvalue != Inf)
},.to="by_site_pvalue") %>>%
  unnest()



##################################  HWE  ################################################
#.row=hwe_test=variatn_lung %>>%filter(ac_lung>2,ac_1kg>2)%>>%{.[10,]}

hwe_test=variatn_lung %>>%
  filter(ac_lung>2,ac_1kg>2) %>>%
  #filter(p_valu<0.0001) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%#head(10)%>>%
  by_row(function(.row){
    .test_tabel=read_tsv("maf_norm/lung/depth/gender_file.tsv") %>>%
      left_join(maf_norm_lung %>>%
                  filter(chr==.row$chr,start==.row$start,end==.row$end,alt==.row$alt)%>>%
                  count(patient_id)) %>>%
      mutate(alt_count=ifelse(is.na(n),0,n)) %>>%
      count(alt_count) %>>%
      dplyr::arrange(alt_count) #これで(refhomo,refalthetero,althomo)の順
    alt_allele_freq=.row$ac_lung / .row$an_lung
    .hwe=c((1-alt_allele_freq) ^2,2*(1-alt_allele_freq) * alt_allele_freq,alt_allele_freq ^2)
    #print(c(.test_tabel$nn,.hwe))
    if(length(.test_tabel$alt_count)==2){#althomoがない場合
      .chisq_test=chisq.test(c(.test_tabel$nn,0),p=.hwe)
      .chisq_test$p.value
    }else{
      .chisq_test=chisq.test(.test_tabel$nn,p=.hwe)
      .chisq_test$p.value
    }
  },.to="hwe_pvalue") %>>%
  unnest()

hwe_test %>>%
  ggplot(aes(x=p_value,y=hwe_pvalue))+geom_point()

############################### PolyPhen 数値 ##################################
variatn_lung %>>%
  filter(ac_lung>2,ac_1kg>2) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  mutate(polyphen_score=str_extract(PolyPhen,"\\d.\\d+")) %>>%
  mutate(polyphen_score=as.double(ifelse(is.na(polyphen_score),str_extract(PolyPhen,"\\d"),polyphen_score))) %>>%
  #mutate(polyphen_score=ifelse((mutype=="splice"|mutype=="truncate"),1,polyphen_score)) %>>%
  filter(!is.na(polyphen_score)) %>>%#View
  ggplot(aes(x=compar,y=polyphen_score))+
  geom_point()


#######################like GWAS plot###############################
if(0){
maf_1kg_all=maf_list %>>%
  mutate(maf=purrr::map(path,~strip_maf(.))) %>>%
  select(-path) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,
                end=End_Position,ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  dplyr::select(sample,allele,gene_symbol,chr,start,end,ref,alt,mutype,PolyPhen,Consequence) %>>%
  left_join(have_y_sample) %>>%
  filter(is.na(Y)) %>>%
  dplyr::select(-Y)
}
#write_df(maf_1kg_all,"/working/1000genomes/maf/extract/all.maf")
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",col_types = "cccciiccccc")

extract_norm_maf_all=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2','Consequence','PolyPhen')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      classify_consequence() %>>%
      filter(Consequence != "intron_variant")
  }
  
  #.bp is bodypart of it cancer
  maf=read_tsv(paste('maf_norm/',.bp,'/depth/gender_file.tsv',sep="")) %>>%
    mutate(filename=paste('maf_norm/',.bp,"/",patient_id,".maf",sep="")) %>>%
    mutate(purrr::map(filename,~strip_maf(.))) %>>%
    unnest() %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,allele1=Tumor_Seq_Allele1,allele2=Tumor_Seq_Allele2) %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand)) %>>%
    tidyr::gather(allele,alt,allele1,allele2) %>>%
    filter(!ref==alt)%>>%
    filter(!(gender=="male"&chr=="X"&allele=="allele2"))%>>%
    dplyr::select(patient_id,gene_symbol,allele,chr,start,end,ref,alt,mutype,PolyPhen,Consequence)
  
}

#lung_norm_maf_all=extract_norm_maf_all("lung")
#write_df(lung_norm_maf_all,"./maf_norm/united_maf/lung_all_site.maf")
lung_norm_maf_all=read_tsv("./maf_norm/united_maf/lung_all_site.maf",col_types = "cccciiccccc")

#breast_norm_maf_all=extract_norm_maf_all("breast")
#write_df(breast_norm_maf_all,"./maf_norm/united_maf/breast_all_site.maf")
breast_norm_maf_all=read_tsv("./maf_norm/united_maf/breast_all_site.maf",col_types = "cccciiccccc")

#brain_norm_maf_all=extract_norm_maf_all("brain")
#write_df(brain_norm_maf_all,"./maf_norm/united_maf/brain_all_site.maf")
brain_norm_maf_all=read_tsv("./maf_norm/united_maf/brain_all_site.maf",col_types = "cccciiccccc")

#kidney_norm_maf_all=extract_norm_maf_all("kidney")
#write_df(kidney_norm_maf_all,"./maf_norm/united_maf/kidney_all_site.maf")
kidney_norm_maf_all=read_tsv("./maf_norm/united_maf/kidney_all_site.maf",col_types = "cccciiccccc")

#colorectal_norm_maf_all=extract_norm_maf_all("colorectal")
#write_df(colorectal_norm_maf_all,"./maf_norm/united_maf/colorectal_all_site.maf")
colorectal_norm_maf_all=read_tsv("./maf_norm/united_maf/colorectal_all_site.maf",col_types = "cccciiccccc")

lung_coverage      =read_tsv("/Volumes/areca42TB/tcga/maf_norm/lung/depth/coverage_all_data_exist_patient.tsv",
                        col_names = c("chr","start","an_lung"),col_types = "cdd")
brain_coverage     =read_tsv("/Volumes/areca42TB/tcga/maf_norm/brain/depth/coverage_all_data_exist_patient.tsv",
                        col_names = c("chr","start","an_brain"),col_types = "cdd")
breast_coverage    =read_tsv("/Volumes/areca42TB/tcga/maf_norm/breast/depth/coverage_all_data_exist_patient.tsv",
                         col_names = c("chr","start","an_breast"),col_types = "cdd")
kidney_coverage    =rad_tsv("/Volumes/areca42TB/tcga/maf_norm/kidney/depth/coverage_all_data_exist_patient.tsv",
                        col_names = c("chr","start","an_breast"),col_types = "cdd")
colorectal_coverage=rad_tsv("/Volumes/areca42TB/tcga/maf_norm/colorectal/depth/coverage_all_data_exist_patient.tsv",
                        col_names = c("chr","start","an_breast"),col_types = "cdd")
cancer_coverage=full_join(lung_coverage,brain_coverage) %>>%
  full_join(breast_coverage) %>>%full_join(kidney_coverage) %>>%
  mutate(an_lung  =ifelse(is.na(an_lung  ),0,an_lung  ),
         an_brain =ifelse(is.na(an_brain ),0,an_brain ),
         an_breast=ifelse(is.na(an_breast),0,an_breast)) %>>%
  mutate(an_cancer=an_lung + an_brain + an_breast) %>>%
  dplyr::select(chr,start,an_cancer)


tally_1kg_all= maf_1kg_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_1kg=n()) %>>%
  left_join(maf_1kg_all %>>% group_by(sample,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_1kg=n())) %>>%
  full_join(vcf_1kg_ac0)
tally_1kg_all=tally_1kg_all %>>%full_join(vcf_1kg_ac0)

tally_lung_all = lung_norm_maf_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_lung=n()) %>>%
  left_join(lung_norm_maf_all %>>% group_by(patient_id,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_lung=n()))
tally_brain_all = brain_norm_maf_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_brain=n()) %>>%
  left_join(brain_norm_maf_all %>>% group_by(patient_id,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_brain=n()))
tally_breast_all = breast_norm_maf_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_breast=n()) %>>%
  left_join(breast_norm_maf_all %>>% group_by(patient_id,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_breast=n()))
tally_kidney_all = breast_kidney_maf_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_breast=n()) %>>%
  left_join(breast_norm_maf_all %>>% group_by(patient_id,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_breast=n()))
tally_colorectal_all = breast_colorectal_maf_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_breast=n()) %>>%
  left_join(breast_norm_maf_all %>>% group_by(patient_id,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_breast=n()))
tally_all=full_join(tally_lung_all,tally_brain_all) %>>%
  full_join(tally_breast_all) %>>%full_join(tally_kidney_all) %>>%full_join(tally_colorectal_all) %>>%
  mutate(ac_lung  =as.double(ifelse(is.na(ac_lung  ),0,ac_lung)  ),
         ac_brain =as.double(ifelse(is.na(ac_brain ),0,ac_brain) ),
         ac_breast=as.double(ifelse(is.na(ac_breast),0,ac_breast))) %>>%
  mutate(homo_lung  =as.double(ifelse(is.na(homo_lung  ),0,homo_lung)  ),
         homo_brain =as.double(ifelse(is.na(homo_brain ),0,homo_brain) ),
         homo_breast=as.double(ifelse(is.na(homo_breast),0,homo_breast))) %>>%
  mutate(ac_cancer= ac_lung + ac_brain + ac_breast ) %>>%
  mutate(homo_cancer= homo_lung + homo_brain + homo_breast) 



test_fisher =function(.row){
  .test=fisher.test(matrix(c(.row$ac_1kg_,.row$ac_cancer,.row$an_1kg - .row$ac_1kg,.row$an_cancer - .row$ac_cancer),2))
  .test$p.value
}

HWE_test=function(AF,homo,hetero,an){
  .hw=c(AF^2,2*AF*(1-AF),1-(AF^2)-(2*AF*(1-AF)))
  .test=chisq.test(c(homo,hetero,(an/2 - homo - hetero)),p=.hw)
  .test$p.value
}
do_HWE=function(.row){
  .cancer=ifelse(is.na(.row$hetero_cancer)|.row$ac_cancer==0,1,
                      HWE_test(.row$AF_cancer,.row$homo_cancer/2,(.row$ac_cancer-.row$homo_cancer),.row$an_cancer))
  .seng=ifelse(is.na(.row$hetero_1kg)|.row$ac_1kg==0,1,
                   HWE_test(.row$AF_1kg,.row$homo_1kg/2,(.row$ac_1kg - .row$homo_1kg),.row$an_1kg))
  .uk10k=ifelse(is.na(.row$hetero_uk10k)|.row$uk_ac==0,1,
                     HWE_test(.row$AF_uk10k,.row$uk_althomo,.row$uk_hetero,.row$uk_an))
  data_frame(HWE_p_cancer=.cancer,HWE_p_1kg=.seng,HWE_p_uk10k=.uk10k)
}

all_variant = full_join(tally_1kg_all,tally_all) %>>%
  mutate(ac_1kg_=as.double(ifelse(is.na(ac_1kg),0,ac_1kg)),
         ac_cancer=as.double(ifelse(is.na(ac_cancer),0,ac_cancer)),
         homo_1kg=as.double(ifelse(is.na(homo_1kg),0,homo_1kg)),
         homo_cancer=as.double(ifelse(is.na(homo_cancer),0,homo_cancer))) %>>%
  left_join(cancer_coverage) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_1kg=ifelse(chr =="X",sample_num_X_1kg,sample_num_1kg)) %>>%
  mutate(hetero_1kg=ifelse(chr=="X",NA,(ac_1kg - homo_1kg)/(an_1kg/2)),
         hetero_cancer=ifelse(chr=="X",NA,(ac_cancer - homo_cancer)/(an_cancer/2))) %>>%
  by_row(.to="p_value",~test_fisher(.)) %>>% unnest() %>>%
  left_join(vcf_10k %>>%filter(nchar(ref)==1,nchar(alt)==1,ref!="-",alt!="-")) %>>%
  left_join(vcf_exac%>>%filter(nchar(ref)==1,nchar(alt)==1,ref!="-",alt!="-") %>>%rename(start=posi)) %>>%
  left_join(vcf_1kg_ac0 %>>%dplyr::select(chr,start,ref,alt,an_1kg)%>>%mutate(focal_1kg_ac0="ac0")) %>>%
  mutate(ac_1kg=ifelse(is.na(focal_1kg_ac0),ac_1kg,0)) %>>%
  mutate(hetero_uk10k =ifelse(chr=="X",NA,ifelse(is.na(uk_an),NA,uk_hetero/(uk_an/2))),
         AF_uk10k=ifelse(is.na(uk_an),NA,uk_ac /uk_an),
         AF_exac=ifelse(is.na(an_exac),NA,ac_exac / an_exac),
         AF_1kg=ac_1kg / an_1kg,
         AF_cancer=ac_cancer/ an_cancer) %>>%
  mutate(
    HWhet_cancer=ifelse(is.na(hetero_cancer),NA,2*AF_cancer*(1 - AF_cancer)),
    HWhet_1kg=ifelse(is.na(hetero_1kg),NA,2*AF_1kg*(1-AF_1kg)),
    HWhet_uk10k=ifelse(is.na(hetero_uk10k),NA,2*AF_uk10k*(1-AF_uk10k)))%>>%
  purrr::by_row(.to="HWE",do_HWE) %>>%unnest()%>>%
  dplyr::select(gene_symbol,chr,start,ref,alt,mutype,ac_cancer,an_cancer,ac_1kg,an_1kg,
                AF_cancer,AF_1kg,AF_uk10k,AF_exac,hetero_cancer,HWhet_cancer,HWE_p_cancer,
                hetero_1kg,HWhet_1kg,HWE_p_1kg,hetero_uk10k,HWhet_uk10k,HWE_p_uk10k,p_value)

all_hwe_variant=all_variant %>>%
  filter(!(HWE_p_1kg<0.01|HWE_p_cancer<0.01|HWE_p_uk10k<0.01))
write_df(all_hwe_variant%>>%dplyr::arrange(p_value),"all_variant.tsv")

all_hwe_variant %>>% filter(p_value >1e-10) %>>% 
  filter(ac_cancer>2,ac_1kg>2) %>>%#filter(gene_symbol=="KMT2C") %>>%
  ggplot()+
  geom_point(aes(x=start, y=p_value, fill= mutype),shape=21)+
  scale_y_log10()+
  scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue",silent="white"))


all_hwe_variant %>>% filter(mutype != "silent") %>>%
  #filter(p_value<0.0001) %>>%
  ggplot()+
  geom_histogram(aes(x=p_value,fill=mutype),binwidth = 0.01)

lung_norm_maf_all %>>%filter(chr==,start==)
##########################################################################
####  homo の数で比べてみたら？　#####
test_fisher_1kg =function(.row){
  .test=fisher.test(matrix(c(.row$homo_1kg/2,.row$homo_cancer/2,
                             (.row$an_1kg - .row$homo_1kg)/2,(.row$an_cancer - .row$homo_cancer)/2),2))
  .test$p.value
}
test_fisher_uk =function(.row){
  .test=fisher.test(matrix(c(.row$uk_althomo,.row$homo_cancer/2,
                             .row$uk_an/2 -.row$uk_althomo,(.row$an_cancer - .row$homo_cancer)/2),2))
  .test$p.value
}

all_variant_homo_fisher = full_join(tally_1kg_all,tally_all) %>>%
  mutate(ac_1kg_=as.double(ifelse(is.na(ac_1kg),0,ac_1kg)),
         ac_cancer=as.double(ifelse(is.na(ac_cancer),0,ac_cancer)),
         homo_1kg=as.double(ifelse(is.na(homo_1kg),0,homo_1kg)),
         homo_cancer=as.double(ifelse(is.na(homo_cancer),0,homo_cancer))) %>>%
  left_join(cancer_coverage) %>>%
  filter(!is.na(an_cancer)) %>>%filter(an_cancer!=0) %>>%
  mutate(an_1kg=ifelse(chr =="X",sample_num_X_1kg,sample_num_1kg)) %>>%
  mutate(hetero_1kg=ifelse(chr=="X",NA,(ac_1kg - homo_1kg)/(an_1kg/2)),
         hetero_cancer=ifelse(chr=="X",NA,(ac_cancer - homo_cancer)/(an_cancer/2))) %>>%
  filter(chr!="X") %>>%
  left_join(vcf_10k %>>%filter(nchar(ref)==1,nchar(alt)==1,ref!="-",alt!="-")) %>>%
  left_join(vcf_exac%>>%filter(nchar(ref)==1,nchar(alt)==1,ref!="-",alt!="-") %>>%rename(start=posi)) %>>%
  left_join(vcf_1kg_ac0 %>>%dplyr::select(chr,start,ref,alt,an_1kg)%>>%mutate(focal_1kg_ac0="ac0")) %>>%
  mutate(uk_althomo=ifelse(is.na(uk_althomo),0,uk_althomo),
         uk_an=ifelse(is.na(uk_an),0,uk_an),
         homo_1kg=ifelse(is.na(homo_1kg),0,homo_1kg)) %>>%
  filter(!(homo_cancer==0 & homo_1kg==0)) %>>% filter(!(homo_cancer==0 & uk_althomo==0)) %>>%
  by_row(.to="p_value_1kg",~test_fisher_1kg(.)) %>>%
  by_row(.to="p_value_uk10k",~test_fisher_uk(.)) %>>% unnest() %>>%
  mutate(ac_1kg=ifelse(is.na(focal_1kg_ac0),ac_1kg,0)) %>>%
  mutate(hetero_uk10k =ifelse(chr=="X",NA,ifelse(is.na(uk_an),NA,uk_hetero/(uk_an/2))),
         AF_uk10k=ifelse(is.na(uk_an),NA,uk_ac /uk_an),
         AF_exac=ifelse(is.na(an_exac),NA,ac_exac / an_exac),
         AF_1kg=ac_1kg / an_1kg,
         AF_cancer=ac_cancer/ an_cancer) %>>%
  mutate(
    HWhet_cancer=ifelse(is.na(hetero_cancer),NA,2*AF_cancer*(1 - AF_cancer)),
    HWhet_1kg=ifelse(is.na(hetero_1kg),NA,2*AF_1kg*(1-AF_1kg)),
    HWhet_uk10k=ifelse(is.na(hetero_uk10k),NA,2*AF_uk10k*(1-AF_uk10k)))%>>%
  purrr::by_row(.to="HWE",do_HWE) %>>%unnest()%>>%
  dplyr::select(gene_symbol,chr,start,ref,alt,mutype,PolyPhen,homo_cancer,homo_1kg,uk_althomo,
                p_value_1kg,p_value_uk10k,AF_cancer,AF_1kg,AF_uk10k,AF_exac,hetero_cancer,
                HWhet_cancer,HWE_p_cancer,hetero_1kg,HWhet_1kg,HWE_p_1kg,hetero_uk10k,HWhet_uk10k,HWE_p_uk10k)

all_variant_homo_fisher %>>%
  filter(mutype=="silent"|mutype=="missense") %>>%
  ggplot()+
  geom_histogram(aes(x=p_value_1kg))+
  facet_grid(~ mutype)
