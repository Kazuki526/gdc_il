library(tidyr)
#library(plyr)
library(dplyr)
library(pipeR)
library(stringr)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(readr)
library(readxl)
library(XML)
library(gtools)
library(purrr)
library(purrrlyr)
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
topdriver_bed=read_tsv("/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed",
                       col_names = c("chr","gene_start","gene_end","ids","score","strand"))%>>%
  tidyr::separate(ids,c("gene_symbol","ensg"),sep=";")
classify_consequence = function(.data) {
  mutate(.data,
         mutype= dplyr::recode(Consequence,
                               downstream_gene_variant = 'flank',
                               `3_prime_UTR_variant` = 'flank',
                               `3_prime_UTR_variant,NMD_transcript_variant`='NMD',
                               upstream_gene_variant = 'flank',
                               `5_prime_UTR_variant` = 'flank',
                               frameshift_variant = 'truncating',
                               frameshift_variant = 'truncating',
                               inframe_deletion = 'inframe_indel',
                               inframe_insertion = 'inframe_indel',
                               intron_variant = 'silent',
                               splice_region_variant = 'splice_region',
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
                               `coding_sequence_variant,3_prime_UTR_variant`='flank',
                               `coding_sequence_variant,5_prime_UTR_variant`='flank',
                               `stop_gained,frameshift_variant` = 'truncating',
                               `stop_gained,frameshift_variant,start_lost` = 'truncating',
                               `missense_variant,splice_region_variant`='missense',
                               `intron_variant,non_coding_transcript_variant`='silent',
                               `intron_variant,NMD_transcript_variant`='NMD',
                               `non_coding_transcript_exon_variant,non_coding_transcript_variant`='silent',
                               `frameshift_variant,splice_region_variant`='truncating',
                               `frameshift_variant,start_lost`='truncating',
                               `frameshift_variant,stop_lost`='truncating',
                               `inframe_deletion,splice_region_variant`='splice_region',
                               `inframe_insertion,splice_region_variant`='splice_region',
                               `frameshift_variant,stop_lost,splice_region_variant`='truncating',
                               `frameshift_variant,stop_retained_variant`='inframe_indel',
                               `protein_altering_variant,splice_region_variant`='inframe_indel',
                               `splice_acceptor_variant,intron_variant`='splice',
                               `splice_acceptor_variant,coding_sequence_variant`='splice',
                               `splice_acceptor_variant,coding_sequence_variant,intron_variant`='splice',
                               `splice_acceptor_variant,5_prime_UTR_variant`='splice',
                               `splice_donor_variant,coding_sequence_variant`='splice',
                               `splice_donor_variant,coding_sequence_variant,intron_variant`='splice',
                               `splice_donor_variant,intron_variant`='splice',
                               `splice_donor_variant,3_prime_UTR_variant,intron_variant`='splice',
                               `splice_donor_variant,5_prime_UTR_variant`='splice',
                               `splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant`='splice',
                               `splice_region_variant,5_prime_UTR_variant`='flank',
                               `splice_region_variant,3_prime_UTR_variant`='flank',
                               `splice_region_variant,intron_variant` = 'splice_region',
                               `splice_region_variant,synonymous_variant`='silent',
                               `splice_region_variant,stop_retained_variant`='silent',
                               `start_lost,inframe_deletion`='truncating',
                               `start_lost,splice_region_variant`='truncating',
                               `stop_lost,inframe_deletion`='truncating',
                               `stop_gained,protein_altering_variant`='truncating',
                               `stop_gained,frameshift_variant,splice_region_variant`='truncating',
                               `stop_gained,inframe_deletion`='truncating',
                               `stop_gained,inframe_insertion`='truncating',
                               `stop_gained,splice_region_variant`='truncating'))
}
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  mutate(role=ifelse(role=="oncogene/TSG","TSG",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))

all_patient_info=read_tsv("~/git/all_patient/all_patient_response.tsv") %>>%
  dplyr::rename(patient_id=submitter_id,age=diagnoses.0.age_at_diagnosis,gender=demographic.gender,
                race=demographic.race,ethnicity=demographic.ethnicity)
patient_race = all_patient_info %>>%
  mutate(race_=ifelse(race=="white" &ethnicity!="hispanic or latino","white",
                      ifelse(race=="black or african american" &ethnicity!="hispanic or latino",
                             "black","other"))) %>>%
  mutate(race_ = ifelse(is.na(race_),"other",race_)) %>>%
  dplyr::select(patient_id,race_) %>>%
  dplyr::rename(race=race_)

########################################################
################ germline mutation #####################
########################################################
if(0){
extract_norm_maf=function(.bp){
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2',"Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2",
                'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position",
                't_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','Transcript_ID')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c', 'c','c','c', 'c','c','c','c','c', 'd','d','d','d','d','d','c'), .)} %>>%
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
    mutate(soma_or_germ=ifelse((t_allele1==n_allele1)&(t_allele2==n_allele2),"germline","somatic"),
           t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
           n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
    mutate(LOH=ifelse((t_genotype=="homo")&(n_genotype=="hetero"),"LOH","no"))
  
}
  data.frame(body_part=c("breast","brain","lung","kidney","colorectal")) %>>%
    mutate(purrr::map(body_part,~extract_norm_maf(.))) %>>%
    unnest() %>>%
    write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/body_part_all.maf.gz")
  
  data.frame(cancer_type=c("hnsc","ov","prad","thca","ucec")) %>>%
    mutate(purrr::map(cancer_type,~extract_norm_maf(.))) %>>%
    unnest()%>>%
    write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/cancer_type_all.maf.gz")
}

#body_part_maf=read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/body_part_all.maf.gz") %>>%
body_part_maf=read_tsv("/Volumes/areca42TB/tcga/all_patient/body_part_all.maf.gz") %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  dplyr::select(-filename) %>>%
  left_join(all_patient_info %>>%dplyr::select(patient_id,cancer_type))

#cancer_type_maf=read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/cancer_type_all.maf.gz") %>>%
cancer_type_maf=read_tsv("/Volumes/areca42TB/tcga/all_patient/cancer_type_all.maf.gz") %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  dplyr::select(-filename) %>>%
  mutate(cancer_type = toupper(cancer_type))

norm_maf_all = rbind(body_part_maf%>>%select(-body_part),cancer_type_maf) %>>%
  mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH)) %>>%
  mutate(cancer_type = ifelse((cancer_type=="COAD" | cancer_type=="READ"),"CRC",cancer_type))%>>%
  mutate(cancer_type = ifelse((cancer_type=="KIRC" | cancer_type=="KIRP" | cancer_type=="KICH"),"KCC",cancer_type))
rm(body_part_maf,cancer_type_maf)

####################################
### read coverage files ####
coverage_all_by_cancer_type =
  read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_all.tsv.gz")
coverage_all = coverage_all_by_cancer_type %>>%
  mutate(an_cancer = an_white + an_black + an_other) %>>%
  group_by(chr,start) %>>%
  summarise(an_cancer = sum(an_cancer))

coverage_male_x_by_cancer_type =
  read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_X_male.tsv.gz") %>>%
  rename(an_white_male = an_white, an_black_male = an_black, an_other_male = an_other)
coverage_male_x = coverage_male_x_by_cancer_type %>>%
  mutate(an_male_cancer = an_white_male + an_black_male + an_other_male) %>>%
  group_by(chr,start) %>>%
  summarise(an_male_cancer = sum(an_male_cancer))

tally_norm_maf = norm_maf_all%>>%
  #filter(cancer_type != "KICH") %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(!(chr=="chrX" & gender=="male" & allele=="n_allele1")) %>>%
  mutate(homo=ifelse(n_genotype=="homo",ifelse(chr=="X" & gender=="male",0,1),0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
            cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
  ungroup() %>>%
  left_join(coverage_all) %>>%left_join(coverage_male_x) %>>%
  mutate(an_male_cancer = ifelse(is.na(an_male_cancer),0,an_male_cancer)) %>>%
  mutate(an_cancer = an_cancer - an_male_cancer) %>>%dplyr::select(-an_male_cancer)


#####################################################################################################################
#### patient number ####
norm_maf_all%>>%
  dplyr::count(cancer_type,patient_id) %>>%
  dplyr::select(-n) %>>%
  dplyr::count(cancer_type) %>>%
  #KICH=66 so delete
  #filter(cancer_type != "KICH") %>>%
  (?sum(.$n))

### patient list ###
patient_list = norm_maf_all %>>%
  dplyr::count(cancer_type,patient_id,age,gender) %>>%
  dplyr::select(-n) %>>%
  filter(!is.na(age))
#write_df(patient_list,"/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
#### AF > 5% のsiteのcoverage(patientごと) #########
## AF_mid_list.tsvはvarscan_maf_age_plot.Rにて作成
mid_af_coverage =read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_coverage_by_patient.tsv.gz") %>>%
  left_join(patient_list) %>>%
  filter(!is.na(cancer_type),focal!="no") %>>%
  dplyr::select(-focal) 



####################################################################################################################
########## somatic mutation ############
#mutation抽出はgerm_pvalu_plot.Rにて
read_maf = function(.file_name){
  read_tsv(.file_name) %>>%
    mutate(patient_id = str_extract(Tumor_Sample_Barcode,"TCGA-\\w\\w-\\w\\w\\w\\w")) %>>%
    dplyr::select(-Tumor_Sample_Barcode,-strand,-cancer_type,-FILTER)
}
#somatic mutationのmaf
somatic_maf = data.frame(cancer_type = c("BRCA","COAD","GBM","HNSC","KICH","KIRC","KIRP","LGG",
                                         "LUAD","LUSC","OV","PRAD","READ","THCA","UCEC")) %>>%
  mutate(file=paste0("kaz_maf/extracted_maf/",tolower(cancer_type),"_topdriver105genes.maf"))　%>>%
  mutate(data = purrr::map(file,~read_maf(.))) %>>%
  unnest() %>>%dplyr::select(-file) %>>%
  mutate(cancer_type = ifelse((cancer_type=="COAD" | cancer_type=="READ"),"CRC",as.character(cancer_type))) %>>%
  #mafはおなじpatientから複数回sequenceしている場合複数回lineがあるため取り除く
  group_by(gene_symbol,chr,start,ref,allele1,allele2,patient_id) %>>%
  filter(t_depth == max(t_depth)) %>>%ungroup() %>>%
  (?.%>>%dplyr::count(gene_symbol,Transcript_ID)%>>%group_by(gene_symbol)%>>%filter(n() >1))
## ENST00000507379とENST00000300305,ENST00000610664は除く
somatic_maf = somatic_maf %>>%
  filter(Transcript_ID !="ENST00000507379",Transcript_ID != "ENST00000300305", Transcript_ID != "ENST00000610664")
# これで１つのgene に1つのtranscript_IDとなる。

###################################################################################################
########## CNA data ##########
#ascat data
extract_ascat = function(.bp){
  read_tsv(paste0("/Volumes/areca42TB/tcga/CNA/",.bp,"/cel/annotate_ascat.tsv.gz"),
           col_types = "ccddcdddddddd")
}
bp_ascat = data.frame(bp = c("breast","brain","lung","kidney","colorectal")) %>>%
  mutate(purrr::map(bp,~extract_ascat(.))) %>>%
  unnest() %>>%
  left_join(all_patient_info %>>%dplyr::select(patient_id,cancer_type)) %>>%
  dplyr::select(-bp)
cantype_ascat = data.frame(cancer_type = c("hnsc","ov","prad","thca","ucec")) %>>%
  mutate(purrr::map(cancer_type, ~extract_ascat(.))) %>>%
  mutate(cancer_type = toupper(cancer_type)) %>>%
  unnest()

annotate_ascat=rbind(bp_ascat,cantype_ascat)%>>%
  mutate(cancer_type = ifelse((cancer_type=="COAD" | cancer_type=="READ"),"CRC",cancer_type))
rm(bp_ascat,cantype_ascat)


#####################################################################################################
########################################### ranking #################################################
#片落ちranking
cna_ranking = annotate_ascat %>>%
  group_by(cancer_type,patient_id,gene_symbol) %>>%
  summarise(nmajor=min(nmajor),nminor=min(nminor)) %>>% ungroup() %>>%
  filter(nminor==0) %>>%
  count(cancer_type,gene_symbol) %>>% ungroup() %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  group_by(cancer_type) %>>%
  mutate(loh_rank = min_rank(desc(n))) 

#somatic snpのranking
somatic_ranking = somatic_maf %>>%
  tidyr::separate(CDS_position,into = c("cds_position","cds_length")) %>>%
  group_by(cancer_type,gene_symbol,mutype) %>>%
  summarise(n=n(),cds_length=as.double(first(cds_length))) %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  group_by(cancer_type,mutype,role) %>>%
  mutate(somatic_rank = min_rank(desc(n/cds_length)))%>>%
  dplyr::select(-n)
