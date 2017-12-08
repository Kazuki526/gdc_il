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

all_patient_info=read_tsv("~/git/all_patient/all_patient_response.tsv") %>>%
  dplyr::rename(patient_id=submitter_id)

########################################################
################ germline mutation #####################
########################################################
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
if(0){
  data.frame(body_part=c("breast","brain","lung","kidney","colorectal")) %>>%
    mutate(purrr::map(body_part,~extract_norm_maf(.))) %>>%
    unnest() %>>%
    write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/body_part_all.maf.gz")
  
  data.frame(cancer_type=c("hnsc","ov","prad","thca","ucec")) %>>%
    mutate(purrr::map(cancer_type,~extract_norm_maf(.))) %>>%
    unnest()%>>%
    write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/cancer_type_all.maf.gz")
}

body_part_maf=read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/body_part_all.maf.gz") %>>%
  filter(mutype!="flank") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  dplyr::select(-filename) %>>%
  left_join(all_patient_info %>>%dplyr::select(patient_id,cancer_type))

cancer_type_maf=read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/cancer_type_all.maf.gz") %>>%
  filter(mutype!="flank") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  dplyr::select(-filename) %>>%
  mutate(cancer_type = toupper(cancer_type))

norm_maf = rbind(body_part_maf%>>%select(-body_part),cancer_type_maf) %>>%
  mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH)) %>>%
  mutate(cancer_type = ifelse((cancer_type=="COAD" | cancer_type=="READ"),"CRC",cancer_type)) 
rm(body_part_maf,cancer_type_maf)

#normalでaltalt, tumorでrefaltとなってる際にnormalでrefのdepth=0のものだけ採用！
#また同じサイトでこのエラーが有意に多い(100patient以上の)siteは解析に使用しないことにした。(3site)
varscan_error = norm_maf %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(an_error = n()*2) %>>%
  ungroup()

norm_maf = norm_maf %>>%
  filter(!(LOH == "back_mutation" & n_ref_count !=0 )) %>>%
  left_join(varscan_error %>>% filter(an_error >=100) %>>%rename(n_allele2=alt)) %>>%
  filter(is.na(an_error)) %>>% dplyr::select(-an_error)%>>%
  filter(chr!="chrX")

tally_norm_maf = norm_maf %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(!(chr=="chrX" & gender=="male" & allele=="n_allele1")) %>>%
  mutate(homo=ifelse(n_genotype=="homo",ifelse(chr=="X" & gender=="male",0,1),0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),homo_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
            cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
  ungroup()
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


############################################################################################################
### read coverage files ####
coverage_all = read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_all.tsv")
coverage_male_x = read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_X_male.tsv")

#### AF > 5% のsiteのcoverage(patientごと) #########
mid_af_coverage =read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_coverage_by_patient.tsv") %>>%
  left_join(varscan_error %>>%filter(an_error > 100) %>>%dplyr::select(chr,start,an_error)) %>>%
  filter(is.na(an_error)) %>>%
  dplyr::select(-an_error) %>>%
  left_join(error_region) %>>%
  filter(!(start >= error_start & start <= error_end)|is.na(error_start)) %>>%
  dplyr::select(-error_start,-error_end,-gene_symbol) %>>%
  left_join(patient_list) %>>%
  filter(!is.na(cancer_type),focal!="no") %>>%
  dplyr::select(-focal)

###################################################################################################################
#prepare  EXAC(nonTCGA) data sets
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",col_types = "cdccdd")
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz",col_types = "cdccdd")
vcf_exac =rbind(vcf_exac,vcf_exac_indel)
rm(vcf_exac_indel)


#######################################################################################################
################################filtering duplicate? position #########################################
#######################################################################################################
### HWE test ###
HWE_test=function(AF,homo,hetero,an){
  .hw=c(AF^2,2*AF*(1-AF),1-(AF^2)-(2*AF*(1-AF)))
  .test=chisq.test(c(homo,hetero,(an/2 - homo - hetero)),p=.hw)
  .test$p.value
}
do_HWE = function(.data){
  HWE_test(.data$ac_cancer/.data$an_cancer,
           .data$homo_cancer,.data$hetero_cancer,.data$an_cancer)
}

### altであるのがheteroのみでduplicateしているっぽい箇所(hetero >= 400 && alt_homo <=10のsite) ###
error_site = tally_norm_maf%>>%
  filter(chr != "chrX") %>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  mutate(start = ifelse(alt == "-",start -1,start)) %>>%
  left_join(coverage_all ) %>>% 
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  mutate(hetero_cancer = ac_cancer - homo_cancer*2) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE=purrr::map(data,~do_HWE(.))) %>>%
  unnest() %>>%
  mutate(ref_homo_cancer = (an_cancer/2 - homo_cancer - hetero_cancer)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,HWE,ac_cancer,an_cancer,ref_homo_cancer,hetero_cancer,homo_cancer) %>>%
  filter(HWE < 0.01) %>>%
  filter(hetero_cancer >=400,homo_cancer <= 10)
write_df(error_site,"/Volumes/areca42TB2/gdc/varscan/all_patient/duplicated_error_site.tsv")

## ちょっとクライテリアをかえてみたら、、って見せるよう
tally_norm_maf%>>%
  filter(chr != "chrX") %>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  mutate(start = ifelse(alt == "-",start -1,start)) %>>%
  left_join(coverage_all ) %>>% 
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  mutate(hetero_cancer = ac_cancer - homo_cancer*2) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE=purrr::map(data,~do_HWE(.))) %>>%
  unnest() %>>%
  mutate(ref_homo_cancer = (an_cancer/2 - homo_cancer - hetero_cancer)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,HWE,ac_cancer,an_cancer,ref_homo_cancer,hetero_cancer,homo_cancer) %>>%
  filter(HWE < 0.01) %>>%
  filter(hetero_cancer <400,hetero_cancer >=200,homo_cancer > 10,homo_cancer < 50) %>>%
  write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/not_duplicated_error_site.tsv")

#error_siteの前後100bpを除く!!(KMT2Cは多すぎるので遺伝子ごと解析から外す)
#もちろんsomatic mutationも!!
error_region = error_site %>>%
  mutate(error_start = start -101, error_end = end +100) %>>%
  group_by(gene_symbol) %>>%
  summarise(chr=first(chr),error_start=min(error_start),error_end=max(error_end)) %>>% ungroup() %>>%
  dplyr::select(gene_symbol,chr,error_start,error_end)

norm_maf = norm_maf %>>%
  filter(gene_symbol != "KMT2C") %>>%
  left_join(error_region) %>>%
  filter(!(start >= error_start & end <= error_end)|is.na(error_start)) %>>%
  dplyr::select(-error_start,-error_end)

somatic_maf = somatic_maf %>>%
  filter(gene_symbol != "KMT2C") %>>%
  left_join(error_region) %>>%
  filter(!(start >= error_start & end <= error_end)|is.na(error_start)) %>>%
  dplyr::select(-error_start,-error_end)

#####################################################################################################################
#### patient number ####
norm_maf%>>%
  dplyr::count(cancer_type,patient_id) %>>%
  dplyr::select(-n) %>>%
  left_join(annotate_ascat %>>%dplyr::count(patient_id) %>>%dplyr::select(-n)%>>% mutate(ascat_focal="ok")) %>>%
  dplyr::count(cancer_type,ascat_focal) %>>%
  #KICH=66 so delete
  filter(cancer_type != "KICH") %>>%
  (?sum(.$n))

### patient list ###
patient_list = norm_maf %>>%
  dplyr::count(cancer_type,patient_id,age) %>>%
  dplyr::select(-n) %>>%
  filter(cancer_type != "KICH") %>>%
  filter(!is.na(age))


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

######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
######## classわけ (MAF=)####
#(~5, 5~1, 1~0.5, 0.5~0.05, 0.05~)
# a   b    c      d        e
classed_site=tally_norm_maf %>>%
  left_join(vcf_exac%>>%mutate(chr=paste0("chr",chr))%>>%dplyr::rename(start=posi))%>>%
  left_join(coverage_all ) %>>%
  mutate(AF = ifelse(is.na(ac_exac),ac_cancer/an_cancer,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))

mid_maf = mid_af_coverage %>>%
  left_join(classed_site) %>>%
  dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac) %>>%
  filter(AF > 0.005) %>>%
  filter(gene_symbol != "KMT2C") %>>%
  left_join(error_region) %>>%
  filter(!(start >= error_start & end <= error_end)|is.na(error_start)) %>>%
  dplyr::select(-error_start,-error_end) %>>%
  left_join(norm_maf %>>%select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                                -Protein_position,-strand,-t_genotype,-n_genotype)) %>>%
  mutate(t_allele1 = ifelse(is.na(t_allele1),ref,t_allele1),
         t_allele2 = ifelse(is.na(t_allele2),ref,t_allele2),
         n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
         n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
         soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
         LOH = ifelse(is.na(LOH),"ref",LOH))

a_maf = mid_maf %>>%filter(MAF>=0.05, mutype!="silent")
b_maf = mid_maf %>>%filter(MAF<0.05,MAF>=0.01, mutype!="silent")
c_maf = mid_maf %>>%filter(MAF<0.01,MAF>=0.005, mutype!="silent")

ref_minor_ef = mid_maf %>>%
  filter(AF>0.995) %>>%
  filter(ref==n_allele1) %>>%
  mutate(soma_or_germ = ifelse(soma_or_germ == "ref","ref_minor",soma_or_germ))
rm(mid_maf)

d_maf = norm_maf %>>% 
  dplyr::select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                                    -Protein_position,-strand,-t_genotype,-n_genotype) %>>%
  left_join(classed_site %>>%dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac)) %>>%
  rbind(ref_minor_ef) %>>%
  filter(MAF<0.005,MAF>=0.0005, mutype!="silent",cancer_type!="KICH")
e_maf = norm_maf %>>%
  dplyr::select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                                    -Protein_position,-strand,-t_genotype,-n_genotype) %>>%
  left_join(classed_site %>>%dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac)) %>>%
  rbind(ref_minor_ef) %>>%
  filter(MAF<0.0005, mutype!="silent",cancer_type!="KICH")

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#ここから解析
#ここから解析
#ここから解析

####################################################################################################################
######### cumulative mutation effect ##########
#まずtruncateを２つ以上もつ患者いる？？
truncating_count = norm_maf %>>%
  filter(!(soma_or_germ =="somatic" & LOH=="no"),gene_symbol!="KMT2C",
         cancer_type!="READ",cancer_type!="KICH") %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating",role=="TSG") %>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count=n())

#truncateの数だけでは？？
.plot = patient_list %>>%
  left_join(truncating_count) %>>%
  mutate(truncating_count=ifelse(is.na(truncating_count),0,truncating_count),
         age=round(age/365.25*100)/100) %>>%
  ggplot(aes(x=as.factor(truncating_count), y=age))+
  geom_violin()+
  geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
  facet_wrap( ~ cancer_type, ncol = 5)+
  ylim(0,90)+
  geom_text(data =truncating_count %>>%count(cancer_type,truncating_count),
            aes(x=as.factor(truncating_count),y=5,label=n),size=3,position="stack")
.plot
ggsave("age_plot/cumulative/allgene_truncating.pdf",.plot,height = 15,width = 15)
#遺伝子絞らないとな、、、、

cumulative_plot = function(.tbl,.class){
##missense の数
  missense_count = .tbl %>>%
    dplyr::select(-alt) %>>%
    left_join(norm_maf %>>%filter(mutype=="missense",!(soma_or_germ =="somatic" & LOH=="no"))) %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role=="TSG") %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #mutate(missense_num=missense_num %/% 3) %>>% #a_mafをやる時用
    mutate(missense_count=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num)) 
#相関直線を
  regression = patient_list %>>%
    left_join(missense_count) %>>%
    mutate(missense_num=ifelse(is.na(missense_num),0,missense_num),
           age = round(age/365.25*100)/100) %>>%
    tidyr::nest(-cancer_type) %>>%
    mutate(data = purrr::map(data,~as.data.frame(as.list(coef(lm(age ~ missense_num,data = .)))))) %>>%
    unnest()
  write_df(regression,paste0("age_plot/cumulative/",.class,"regression.tsv"))
##バイオリンプロットで見やすく
  .plot = patient_list %>>%
    left_join(missense_count) %>>%
    mutate(missense_count = ifelse(is.na(missense_count),"0",missense_count),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           age=round(age/365.25*100)/100) %>>%
    left_join(truncating_count) %>>%
    mutate(truncating=ifelse(is.na(truncating_count),"not_have_truncating","have_truncating")) %>>%
    filter(truncating=="not_have_truncating") %>>%
    ggplot(aes(x=reorder(missense_count,missense_count_order), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    facet_wrap( ~ cancer_type, ncol = 5)+
    geom_abline(data = regression %>>%filter(missense_num > 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "blue")+
    geom_abline(data = regression %>>%filter(missense_num <= 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "red")+
    ylim(0,90)+
    geom_text(data =missense_count %>>%count(cancer_type,missense_count),
              aes(x=missense_count,y=5,label=n),size=3,position="stack")+
    theme(strip.text = element_text(size=6),axis.text.x = element_text(angle = -45, hjust = 0))
  .plot
  ggsave(paste0("age_plot/cumulative/",.class,".pdf"),.plot,height = 15,width = 15)
}

mid_count = function(.tbl){
  .tbl %>>%
    mutate(alt_count = ifelse(n_allele2==ref,0,ifelse(n_allele1==ref,1,2))) %>>%
    mutate(MAC =ifelse(AF !=MAF,2-alt_count,alt_count)) %>>%
    dplyr::select(-alt_count)
}

a_maf %>>%mid_count%>>%cumulative_plot("a") #(この場合はcount数の表示をいじる必要がある)
b_maf %>>%mid_count%>>%cumulative_plot("b")
c_maf %>>%mid_count%>>%cumulative_plot("c")
d_maf %>>%mutate(MAC=1) %>>%cumulative_plot("d")
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot("0.05")

#bc
rbind(b_maf,c_maf)%>>%mid_count()%>>%cumulative_plot("bc")
#cd
c_maf %>>%mid_count()%>>%rbind(d_mid %>>%mutate(MAC=1))%>>%cumulative_plot("cd")
#de
rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%cumulative_plot("de")

#bcd
rbind(b_maf,c_maf)%>>%mid_count()%>>%rbind(d_mid %>>%mutate(MAC=1))%>>%cumulative_plot("bcd")
#cde
rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%rbind(c_maf%>>%mid_count())%>>%cumulative_plot("cde")

#bcde
rbind(b_maf,c_maf)%>>%mid_count()%>>%
  rbind(rbind(d_maf,e_maf)%>>%mutate(MAC=1))%>>%cumulative_plot("bcde")
