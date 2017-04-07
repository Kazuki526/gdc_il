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
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  mutate(role=ifelse(role=="oncogene/TSG","oncogene",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))

### all norm maf###
lung_norm_maf_all=read_tsv("./maf_norm/united_maf/lung_all_site.maf",col_types = "cccciiccccc")
breast_norm_maf_all=read_tsv("./maf_norm/united_maf/breast_all_site.maf",col_types = "cccciiccccc")
brain_norm_maf_all=read_tsv("./maf_norm/united_maf/brain_all_site.maf",col_types = "cccciiccccc")

###tally  maf###
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
tally_all=full_join(tally_lung_all,tally_brain_all) %>>%
  full_join(tally_breast_all) %>>%
  mutate(ac_lung  =as.double(ifelse(is.na(ac_lung  ),0,ac_lung)  ),
         ac_brain =as.double(ifelse(is.na(ac_brain ),0,ac_brain) ),
         ac_breast=as.double(ifelse(is.na(ac_breast),0,ac_breast))) %>>%
  mutate(homo_lung  =as.double(ifelse(is.na(homo_lung  ),0,homo_lung)  ),
         homo_brain =as.double(ifelse(is.na(homo_brain ),0,homo_brain) ),
         homo_breast=as.double(ifelse(is.na(homo_breast),0,homo_breast))) %>>%
  mutate(ac_cancer= ac_lung + ac_brain + ac_breast ) %>>%
  mutate(homo_cancer= homo_lung + homo_brain + homo_breast)

###coverage file###
lung_coverage  =read_tsv("/Volumes/areca42TB/tcga/maf_norm/lung/depth/coverage_all_data_exist_patient.tsv",
                         col_names = c("chr","start","an_lung"),col_types = "cdd")
brain_coverage =read_tsv("/Volumes/areca42TB/tcga/maf_norm/brain/depth/coverage_all_data_exist_patient.tsv",
                         col_names = c("chr","start","an_brain"),col_types = "cdd")
breast_coverage=read_tsv("/Volumes/areca42TB/tcga/maf_norm/breast/depth/coverage_all_data_exist_patient.tsv",
                         col_names = c("chr","start","an_breast"),col_types = "cdd")
cancer_coverage=full_join(lung_coverage,brain_coverage) %>>%
  full_join(breast_coverage) %>>%
  mutate(an_lung  =ifelse(is.na(an_lung  ),0,an_lung  ),
         an_brain =ifelse(is.na(an_brain ),0,an_brain ),
         an_breast=ifelse(is.na(an_breast),0,an_breast)) %>>%
  mutate(an_cancer=an_lung + an_brain + an_breast) %>>%
  dplyr::select(chr,start,an_cancer)

### somatic mutation####
lusc_maf=read_tsv("kaz_maf/extracted_maf/lusc_topdriver105genes.maf")
luad_maf=read_tsv("kaz_maf/extracted_maf/luad_topdriver105genes.maf")
lung_maf=full_join(lusc_maf,luad_maf) %>>%
  left_join(read_tsv("maf_norm/lung/depth/gender_file.tsv")) %>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)
breast_maf=read_tsv("kaz_maf/extracted_maf/brca_topdriver105genes.maf")
lgg_maf=read_tsv("kaz_maf/extracted_maf/lgg_topdriver105genes.maf")
gbm_maf=read_tsv("kaz_maf/extracted_maf/gbm_topdriver105genes.maf")
brain_maf=full_join(lgg_maf,gbm_maf) %>>%
  left_join(read_tsv("maf_norm/brain/depth/gender_file.tsv")) %>>%
  filter(!is.na(gender))%>>%
  dplyr::select(-Tumor_Sample_Barcode,-gender)

maf_somatic=full_join(lung_maf,breast_maf) %>>%full_join(brain_maf)

### CNA ###
lung_cna=read_tsv("./CNA/lung/cel/annotate_ascat.tsv.gz")
breast_cna=read_tsv("./CNA/breast/cel/annotate_ascat.tsv.gz")
brain_cna=read_tsv("./CNA/brain/cel/annotate_ascat.tsv.gz")

###1kg###
sample_num_1kg=5008
sample_num_X_1kg=sample_num_1kg - 1233
vcf_1kg_ac0=read_tsv("/working/1000genomes/maf/extract/ac0_variation.tsv")
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",col_types = "cccciiccccc")

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

####UK 10K ###
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_likevcf.tsv.gz",col_types = "cdccdddd")

### EXAC nonTCGA ###
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",col_types = "cdccdd")

######################################################
####################### test #########################
######################################################

### alt homo frequency###
#AFでは全然だめだったので
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

test_fisher =function(.row){
  .test=fisher.test(matrix(c(.row$ac_1kg_,.row$ac_cancer,.row$an_1kg - .row$ac_1kg_,.row$an_cancer - .row$ac_cancer),2))
  .test$p.value
}

all_variant = full_join(tally_1kg_all,tally_all) %>>%
  mutate(ac_1kg_=as.double(ifelse(is.na(ac_1kg),0,ac_1kg)),
         ac_cancer=as.double(ifelse(is.na(ac_cancer),0,ac_cancer)),
         homo_1kg=as.double(ifelse(is.na(homo_1kg),0,homo_1kg)),
         homo_cancer=as.double(ifelse(is.na(homo_cancer),0,homo_cancer))) %>>%
  left_join(cancer_coverage) %>>%
  filter(!is.na(an_cancer),!is.na(gene_symbol)) %>>% filter(an_cancer!=0) %>>%
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
  mutate(HWhet_cancer=ifelse(is.na(hetero_cancer),NA,2*AF_cancer*(1 - AF_cancer)),
         HWhet_1kg=ifelse(is.na(hetero_1kg),NA,2*AF_1kg*(1-AF_1kg)),
         HWhet_uk10k=ifelse(is.na(hetero_uk10k),NA,2*AF_uk10k*(1-AF_uk10k)))%>>%
  purrr::by_row(.to="HWE",do_HWE) %>>%unnest()%>>%
  dplyr::select(gene_symbol,chr,start,ref,alt,mutype,PolyPhen,
                ac_cancer,an_cancer,ac_1kg,an_1kg,uk_ac,uk_an,
                AF_cancer,AF_1kg,AF_uk10k,AF_exac,hetero_cancer,
                HWhet_cancer,HWE_p_cancer,hetero_1kg,HWhet_1kg,
                HWE_p_1kg,hetero_uk10k,HWhet_uk10k,HWE_p_uk10k,p_value)

all_variant=all_variant %>>%
  filter(gene_symbol !="KMT2C")%>>%
  filter(!(HWhet_cancer < hetero_cancer & HWE_p_cancer<1.0e-10)) %>>%
  filter(!(HWhet_1kg < hetero_1kg & HWE_p_1kg<1.0e-10)) %>>%
  filter(!(HWhet_uk10k < hetero_uk10k & HWE_p_uk10k<1.0e-10))
all_variant %>>%filter(uk_ac>5)%>>%count(mutype)


all_variant %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="missense") %>>%
  mutate(polyphen_level=str_extract(PolyPhen,"^\\w+")) %>>%
  mutate(ac_1kg=ifelse(is.na(ac_1kg),0,ac_1kg)) %>>%
  mutate(uk_ac=ifelse(is.na(uk_ac),0,uk_ac),uk_an=ifelse(is.na(uk_an),4090,uk_an)) %>>%
  filter(ac_cancer/an_cancer<0.1,uk_ac/uk_an<0.1) %>>%
  ggplot()+
  geom_abline(intercept = 0,slope = 1)+
  geom_point(aes(x=ac_cancer/an_cancer,y=uk_ac/uk_an,fill=polyphen_level),shape=21)+
  facet_wrap(~ role)+
  scale_y_log10()+
  scale_x_log10()

all_variant %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating") %>>%
  mutate(ac_1kg=ifelse(is.na(ac_1kg),0,ac_1kg)) %>>%
  mutate(uk_ac=ifelse(is.na(uk_ac),0,uk_ac),uk_an=ifelse(is.na(uk_an),4090,uk_an)) %>>%
  #filter(ac_cancer/an_cancer<0.1,uk_ac/uk_an<0.1) %>>%
  ggplot()+
  geom_abline(intercept = 0,slope = 1)+
  geom_point(aes(x=ac_cancer/an_cancer,y=uk_ac/uk_an),shape=21)+
  facet_wrap(~ role)+
  scale_y_log10()+
  scale_x_log10()


all_variant %>>%
  filter(uk_ac >5) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="missense"|mutype=="silent") %>>%
  #mutate(polyphen_level=str_extract(PolyPhen,"^\\w+")) %>>%
  mutate(ac_1kg=ifelse(is.na(ac_1kg),0,ac_1kg)) %>>%
  mutate(uk_ac=ifelse(is.na(uk_ac),0,uk_ac),uk_an=ifelse(is.na(uk_an),4090,uk_an)) %>>%
  #filter(ac_cancer/an_cancer<0.1,uk_ac/uk_an<0.1) %>>%
  ggplot(aes(x=ac_cancer/an_cancer,y=uk_ac/uk_an))+
  geom_abline(intercept = 0,slope = 1)+
  geom_point(aes(fill=role),shape=21)+
  facet_wrap(~ mutype)+
  stat_smooth(method = "lm", se = FALSE)+
  scale_y_log10()+
  scale_x_log10()


#########################################################################################
 ################ germline mutation はがん化に寄与する？##########

#### AF>10%のsite抽出
AF_more10 = all_variant %>>%
  filter(AF_cancer > 0.1) %>>%
  filter(mutype=="missense") %>>%
  filter(an_cancer > 5000)

maf_norm_breast=read_tsv("maf_norm/united_maf/breast_nonsilent.maf")
maf_norm_lung=read_tsv("maf_norm/united_maf/lung_nonsilent.maf")
maf_norm_brain=read_tsv("maf_norm/united_maf/brain_nonsilent.maf")

maf_AF_more10 =full_join(maf_norm_breast,maf_norm_lung) %>>%
  full_join(maf_norm_brain) %>>%
  left_join(AF_more10 %>>%dplyr::select(gene_symbol,chr,start,alt)%>>%mutate(focal="ok")) %>>%
  filter(focal=="ok") %>>% dplyr::select(-focal) %>>%
  group_by(gene_symbol,chr,start,patient_id,ref,alt,PolyPhen) %>>%
  summarise(genotype=ifelse(n()==2,"alt_homo","hetero"))

maf_somatic_for_somatic_test=maf_somatic %>>%dplyr::select(gene_symbol,patient_id,cancer_type,mutype,PolyPhen) %>>%
  mutate(mutation_score=ifelse(mutype=="truncating",2,
                               ifelse(mutype=="splice",1.5,str_extract(PolyPhen,"\\d.\\d+"))))%>>%
  mutate(mutation_score=as.double(ifelse(is.na(mutation_score),str_extract(PolyPhen,"\\d"),
                                         mutation_score)))%>>%
  group_by(gene_symbol,patient_id)%>>%
  filter(mutation_score==max(mutation_score)) %>>%ungroup() %>>%
  rename(mutype_s=mutype) %>>%
  dplyr::select(gene_symbol,patient_id,mutype_s,mutation_score)

full_join(lung_cna,breast_cna) %>>%
  full_join(brain_cna) %>>%
  dplyr::select(gene_symbol,patient_id,nmajor,nminor) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  left_join(AF_more10) %>>%
  filter(!is.na(ref)) %>>%
  dplyr::select(gene_symbol,patient_id,nmajor,nminor,chr,start,ref,alt,AF_cancer,AF_1kg,AF_uk10k,AF_exac) %>>%
  left_join(maf_AF_more10) %>>%
  mutate(genotype=ifelse(is.na(genotype),"ref_homo",genotype)) %>>%
  left_join(maf_somatic_for_somatic_test) %>>%
  mutate(cna_del=ifelse(nminor==0,ifelse(nmajor==0,"homodel","del"),"no"),
         mutype_s=ifelse(is.na(mutype_s),"no",ifelse(mutation_score>0.5,"truncate","benign"))) %>>%
  mutate(somatic_mutation=ifelse(cna_del=="homodel"|(cna_del=="del" & mutype_s=="truncate"),"double_hit",
                                 ifelse(cna_del=="del"|mutype_s=="truncate","one_hit","none"))) %>>%
  #group_by(gene_symbol,chr,start,ref,alt,AF_cancer,AF_1kg,AF_uk10k,AF_exac,genotype,somatic_mutation,PolyPhen)  %>>%
  #mutate(n=n()) %>>%ungroup()%>>%
  mutate(genotype_order=ifelse(genotype=="ref_homo",1,ifelse(genotype=="hetero",2,3))) %>>%
  mutate(somatic_order=ifelse(somatic_mutation=="double_hit",1,ifelse(somatic_mutation=="one_hit",2,3))) %>>%
  mutate(site=as.character(paste(gene_symbol,paste0("chr",chr),start,
                                 paste0(ref,"to",alt),paste0("AF=",floor(AF_cancer*100)),sep=":"))) %>>%
  ggplot()+
  geom_bar(aes(x=reorder(genotype,genotype_order),fill=reorder(as.factor(somatic_mutation),somatic_order)),
           position = "fill")+
  #geom_text(aes(x=genotype,y=,label=n),size=3,hjust=0.5,vjust=3,position="stack")+
  facet_wrap(~ site,ncol=5)+
  theme(strip.text = element_text(size=8))

#### AF<10% の site ###
maf_AF_less10 =full_join(maf_norm_breast,maf_norm_lung) %>>%
  full_join(maf_norm_brain) %>>%
  filter(gene_symbol!="KMT2C") %>>%
  left_join(AF_more10 %>>%dplyr::select(gene_symbol,chr,start,alt)%>>%mutate(focal="ok")) %>>%
  filter(is.na(focal)) %>>% dplyr::select(-focal) %>>%
  group_by(gene_symbol,chr,start,patient_id,ref,alt,mutype,PolyPhen) %>>%
  summarise(genotype=ifelse(n()==2,"alt_homo","hetero")) %>>%
  left_join(all_variant %>>%select(gene_symbol,chr,start,ref,alt,AF_cancer)) %>>%
  filter(AF_cancer<0.1|is.na(AF_cancer)) %>>%
  ungroup() %>>%
  filter(mutype!="splice") %>>%
  mutate(mutation_score=as.numeric(ifelse(mutype=="truncating",2,
                               ifelse(mutype=="inframe_indel",1.5,str_extract(PolyPhen,"\\d.\\d+")))))%>>%
  mutate(mutation_score=as.double(ifelse(is.na(mutation_score),str_extract(PolyPhen,"\\d"),
                                         mutation_score))) %>>%
  mutate(mutation_type=ifelse(mutype=="missense",str_extract(PolyPhen,"^\\w+"),mutype)) %>>%
  filter(!is.na(mutation_type)) %>>%filter(mutation_type!="unknown")
ggplot(maf_AF_less10)+geom_histogram(aes(x=mutation_score))
maf_AF_less10%>>%dplyr::count(mutype)
maf_AF_less10 %>>%group_by(gene_symbol,patient_id) %>>%filter(mutation_score==max(mutation_score)) %>>%
  tally() %>>%filter(n>1) %>>%View ###ある。。。

mutype_order = function(.data) {
  dplyr::mutate(.data, mutation_type_order= dplyr::recode(mutation_type,
                                             "truncating"=1,
                                             "inframe_indel"=2,
                                             "probably_damaging"=3,
                                             "possibly_damaging"=4,
                                             "benign"=5,
                                             "none"=6))
}

full_join(lung_cna,breast_cna) %>>%
  full_join(brain_cna) %>>%
  dplyr::select(gene_symbol,patient_id,nmajor,nminor) %>>%
  group_by(gene_symbol,patient_id) %>>%
  summarise_each(funs(min(.))) %>>%
  #filter(gene_symbol==c("TP53","APC","BRCA1","BRCA2","PIK3R1","CDKN2A","PTEN","CDKN1B","ARID2","CDH1"))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  left_join(maf_AF_less10 %>>%group_by(gene_symbol,patient_id) %>>% #filter(is.na(AF_cancer)|AF_cancer<0.01)%>>%
              filter(mutation_score==max(mutation_score))%>>%summarise_each(funs(.[1]))) %>>%
  left_join(maf_somatic_for_somatic_test) %>>%
  mutate(cna_del=ifelse(nminor==0,ifelse(nmajor==0,"homodel","del"),"no"),
         mutype_s=ifelse(is.na(mutype_s),"no",ifelse(mutation_score>0.5,"truncate","benign"))) %>>%
  mutate(somatic_mutation=ifelse(cna_del=="homodel"|(cna_del=="del" & mutype_s=="truncate"),"double_hit",
                                 ifelse(cna_del=="del"|mutype_s=="truncate","one_hit","none"))) %>>%
  mutate(mutation_type=ifelse(is.na(mutation_type),"none",mutation_type)) %>>%
  mutate(somatic_order=ifelse(somatic_mutation=="double_hit",1,ifelse(somatic_mutation=="one_hit",2,3))) %>>%
  mutype_order() %>>%#View()
  ggplot()+
  geom_bar(aes(x=reorder(mutation_type,mutation_type_order),fill=reorder(as.factor(somatic_mutation),somatic_order)),
           position = "fill")




  