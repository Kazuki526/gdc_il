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
  mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH))
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
  filter(is.na(an_error)) %>>% dplyr::select(-an_error)
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

annotate_ascat=rbind(bp_ascat,cantype_ascat)
rm(bp_ascat,cantype_ascat)

########### CNAのdel,amp頻度のランキング(cancer type ごとに) ###########
#片落ちranking
loh_ranking = annotate_ascat %>>%
  group_by(cancer_type,patient_id,gene_symbol) %>>%
  summarise(nmajor=min(nmajor),nminor=min(nminor)) %>>% ungroup() %>>%
  filter(nminor==0) %>>%
  count(cancer_type,gene_symbol) %>>% ungroup() %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  group_by(cancer_type) %>>%
  mutate(loh_rank = min_rank(desc(n))) 
# homo deletion ranking
homodel_ranking = annotate_ascat %>>%
  group_by(cancer_type,patient_id,gene_symbol) %>>%
  summarise(nmajor=min(nmajor),nminor=min(nminor)) %>>% ungroup() %>>%
  filter(nmajor==0,nminor==0) %>>%
  count(cancer_type,gene_symbol) %>>% ungroup() %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  group_by(cancer_type) %>>%
  mutate(homodel_rank = min_rank(desc(n)))
# 1 copy になっている数ランキング
del_ranking = annotate_ascat %>>%
  group_by(cancer_type,patient_id,gene_symbol) %>>%
  summarise(nmajor=min(nmajor),nminor=min(nminor)) %>>% ungroup() %>>%
  filter(nmajor <=1,nminor==0) %>>%
  count(cancer_type,gene_symbol) %>>% ungroup() %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  filter(role=="TSG") %>>%
  group_by(cancer_type) %>>%
  mutate(del_rank = min_rank(desc(n)))

cna_del_rank = loh_ranking %>>%dplyr::select(-n) %>>%
  left_join(homodel_ranking %>>%dplyr::select(-n)) %>>%
  left_join(del_ranking %>>%dplyr::select(-n))
rm(loh_ranking,homodel_ranking,del_ranking)
#homo delは関係なさそう、、、



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
  (?.%>>%dplyr::count(gene_symbol,Transcript_ID)%>>%group_by(gene_symbol)%>>%filter(n() >1))
## ENST00000507379とENST00000300305,ENST00000610664は除く
somatic_maf = somatic_maf %>>%
  filter(Transcript_ID !="ENST00000507379",Transcript_ID != "ENST00000300305", Transcript_ID != "ENST00000610664")
# これで１つのgene に1つのtranscript_IDとなる。


### somatic mutation頻度の高い遺伝子(cancer typeごとに)
somatic_ranking = somatic_maf %>>%
  tidyr::separate(CDS_position,into = c("cds_position","cds_length")) %>>%
  group_by(cancer_type,gene_symbol,mutype) %>>%
  summarise(n=n(),cds_length=as.double(first(cds_length))) %>>%
  left_join(driver_genes %>>%select(gene,role),by = c("gene_symbol"="gene")) %>>%
  group_by(cancer_type,mutype,role) %>>%
  mutate(somatic_rank = min_rank(desc(n/cds_length)))%>>%
  dplyr::select(-n)



############################################################################################################
### read coverage files ####
coverage_all = read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_all.tsv")
coverage_male_x = read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/coverage_X_male.tsv")



###################################################################################################################
#prepare 1000genomes UK10K EXAC(nonTCGA) data sets
###1kg###
if(0){
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

rm(vcf_1kg_ac0,maf_1kg_all)
####UK 10K ###
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_likevcf.tsv.gz",col_types = "cdccdddd")
}
### EXAC nonTCGA ###
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",col_types = "cdccdd")
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz",col_types = "cdccdd")
vcf_exac =rbind(vcf_exac,vcf_exac_indel)
rm(vcf_exac_indel)


######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
tally_norm_maf = norm_maf %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(!(chr=="chrX" & gender=="male" & allele=="n_allele1")) %>>%
  mutate(homo=ifelse(n_genotype=="homo",ifelse(chr=="X" & gender=="male",0,1),0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),homo_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype)) %>>%
  ungroup()

#### classify each variant AF = ~0.5%(3),0.5%~5%(2),5%~(3) in AF_class colum
classed_site = tally_norm_maf %>>%
  left_join(vcf_exac%>>%mutate(chr=paste0("chr",chr))%>>%dplyr::rename(start=posi)) %>>%
  mutate(AF_class = ifelse(ac_exac/an_exac > 0.05,1,ifelse(ac_exac/an_exac >0.005,2,3)))%>>%
  mutate(AF_class = ifelse(!is.na(AF_class),AF_class,
                           ifelse(ac_cancer >700,1,ifelse(ac_cancer>100,2,3)))) %>>%
  filter(gene_symbol!="KMT2C") %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  ungroup()%>>%(?.%>>%dplyr::count(AF_class))
if(0){
##### save AF_mid_coverage_by_patient.tsv to conduct perl(make_coverage.pl) ####
  classed_site %>>%
    filter(AF_class==1 | AF_class==2) %>>%
    write_df("/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_list.tsv")
}

#######################################################################################################
################################filtering duplicate? position #########################################
#######################################################################################################
###目で見た結果　明らかにおかしくて除くべきsiteはなかった。
if(0){
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


##### chrX以外 #####
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
  filter(homo_cancer > 50)

##### chr X ######
x_error_site = norm_maf %>>%
  filter(chr =="chrX") %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(gender != "male") %>>%
  mutate(homo=ifelse(n_genotype=="homo",1,0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),homo_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype)) %>>%
  ungroup() %>>%
  left_join(coverage_all) %>>%
  left_join(coverage_male_x ) %>>%
  mutate(an_cancer = an_cancer - an_male_cancer) %>>%
  mutate(hetero_cancer = ac_cancer - homo_cancer*2) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE=purrr::map(data,~do_HWE(.))) %>>%
  unnest() %>>%
  mutate(ref_homo_cancer = (an_cancer/2 - homo_cancer - hetero_cancer)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,HWE,ac_cancer,an_cancer,ref_homo_cancer,hetero_cancer,homo_cancer) %>>%
  filter(HWE < 0.01) 

tally_norm_maf %>>%
  left_join(coverage_tbl %>>% mutate(chr = paste0("chr",chr))) %>>%View
  filter(an_cancer >9000)%>>%
  ggplot()+
  geom_histogram(aes(x=an_cancer),binwidth = 100 )

error_site %>>%
  filter( HWE <1E-100) %>>%View
}
#特になさそう、、、

#####################################################################################################################
#### patient number ####
norm_maf%>>%
  dplyr::count(cancer_type,patient_id) %>>%
  dplyr::select(-n) %>>%
  left_join(annotate_ascat %>>%dplyr::count(patient_id) %>>%dplyr::select(-n)%>>% mutate(ascat_focal="ok")) %>>%
  dplyr::count(cancer_type,ascat_focal) %>>%
#KICH=66, READ=151 so delete
  filter(cancer_type != "KICH" & cancer_type != "READ" ) %>>%
  (?sum(.$n))

### patient list ###
patient_list = norm_maf %>>%
  dplyr::count(cancer_type,patient_id,age) %>>%
  dplyr::select(-n) %>>%
  filter(cancer_type != "KICH" & cancer_type != "READ" ) %>>%
  filter(!is.na(age))

#### class 1,2のsiteのcoverage(patientごと) #########
mid_af_coverage =read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/AF_mid_coverage_by_patient.tsv") %>>%
  left_join(varscan_error %>>%filter(an_error > 100) %>>%dplyr::select(chr,start,an_error)) %>>%
  filter(!is.na(an_error)) %>>%
  dplyr::select(-an_error)



##########################################
########## AF >0.5% site(class 1 or 2) ##########
##########################################
class1_plot = function(.data,.cancer_type){
  .mean_age = .data %>>%
    group_by(site,genotype) %>>%
    summarise(mean = round(mean(age)*100)/100) 
  .patient_num = .data %>>%
    dplyr::count(site,genotype)
  .data %>>%
    ggplot(aes(x=reorder(genotype,genotype_order),y=age))+
    geom_violin()+
    geom_boxplot(width=.2,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=1)+
    facet_wrap(~ site,ncol=6)+
    ylim(0,120)+
    geom_signif(comparisons = list(c("refref","refalt")),
                test = "wilcox.test",map_signif_level = T, col="red",y_position = 95)+
    geom_signif(comparisons = list(c("refalt","altalt")),
                test = "wilcox.test",map_signif_level = T, col="red",y_position = 100)+
    theme(strip.text = element_text(size=6),axis.text.x = element_text(angle = -45, hjust = 0))+
    ggtitle(.cancer_type)+
    geom_text(data =.patient_num,aes(x=genotype,y=5,label=n),size=3,position="stack")+
    geom_text(data = .mean_age,aes(x=genotype,y=120,label=mean),size=2,position = "stack",color="red")
}

#####################################################
##### siteごとにgenotypeで発症年齢の上下を見る ######
#####################################################
## all patient ###(AF_class == 1か2にしてclassを変える)
.plot = classed_site %>>%
  filter(mutype!="silent",AF_class == 2,role == "TSG") %>>%
  mutate(patient_focal="focal") %>>%
  left_join(patient_list %>>%mutate(patient_focal="focal")) %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-ac_exac,an_exac,-AF_class,-patient_focal) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no"))) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt")),
         PolyPhen=ifelse(is.na(PolyPhen),"none",str_extract(PolyPhen,"\\w+")),
         age=round(age/365.25*100)/100) %>>%
  mutate(site=paste(gene_symbol,paste(ref,alt,sep=">"),mutype,PolyPhen,sep = ":"),
         genotype_order = ifelse(genotype=="refref",1,ifelse(genotype=="refalt",2,3))) %>>%
  class1_plot("all")
ggsave("age_plot/class2_all.pdf",.plot,width = 10,height = 15)

## plot by cancer_type ##(同じくAF_class == 1か2に変える)
class1_bycancertype_plot = function(.data){
  .CT=first(.data$CT)
  class1_plot(.data,.CT)
}
.plot_byCT = classed_site %>>%
  filter(mutype!="silent",AF_class == 2,role == "TSG") %>>%
  mutate(patient_focal="focal") %>>%
  left_join(patient_list %>>%mutate(patient_focal="focal")) %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-ac_exac,an_exac,-AF_class,-patient_focal) %>>%
  left_join(norm_maf%>>%filter(!(soma_or_germ=="somatic" & LOH =="no"))) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt")),
         PolyPhen=ifelse(is.na(PolyPhen),"none",str_extract(PolyPhen,"\\w+")),
         age=round(age/365.25*100)/100) %>>%
  mutate(site=paste(gene_symbol,paste(ref,alt,sep=">"),mutype,PolyPhen,sep = ":"),
         genotype_order = ifelse(genotype=="refref",1,ifelse(genotype=="refalt",2,3)),
         CT=cancer_type ) %>>%
  nest(-cancer_type) %>>%
  mutate(plot=purrr::map(data,~class1_bycancertype_plot(.)))
ggsave("age_plot/class2_by_cancer_type.pdf",
       gridExtra::marrangeGrob(.plot_byCT$plot,nrow = 1,ncol = 1,top = NULL),
       height = 15, width = 10)
rm(.plot_byCT)



##########################################
######### AF <5% site(class 2,3) #########
##########################################
class23_plot = function(.site_tbl, .title){
  mutype_order = function(.data) {
    dplyr::mutate(.data, mutation_type_order= dplyr::recode(mutation_type,
                                                            "truncating"=1,
                                                            "splice"=2,
                                                            "damaging"=3,
                                                            "probably_damaging"=4,
                                                            "possibly_damaging"=5,
                                                            "benign"=6,
                                                            "none"=7))
  }
  nofilter = norm_maf %>>%
    left_join(.site_tbl %>>%dplyr::rename(n_allele2 = alt)) %>>%
    filter(!is.na(AF_class),!(soma_or_germ=="somatic" & LOH=="no")) %>>%
    mutate(PolyPhen=str_extract(PolyPhen,"^\\w+")) %>>%
    mutate(mutation_type = ifelse(mutype=="truncating",mutype,PolyPhen))%>>%
    mutate(mutation_type = ifelse(str_detect(mutation_type,"damaging"),"damaging",mutation_type)) %>>%
    mutate(mutation_type = ifelse(mutation_type=="unknown","benign",mutation_type),
           site_focal="have") %>>%#(?.%>>%count(mutation_type))%>>%
    mutype_order() %>>%
    dplyr::rename(alt=n_allele2) %>>%
    group_by(patient_id,gene_symbol) %>>%
    filter(mutation_type_order == min(mutation_type_order)) %>>%
    summarise_all(first) %>>%
    dplyr::select(patient_id,gene_symbol,mutation_type,mutation_type_order)

  class23_ggplot = function(.data) {
    .mean_age = .data %>>%
      group_by(cancer_type,mutation_type) %>>%
      summarise(mean = round(mean(age)*100)/100) 
    .patient_num = .data %>>%
      dplyr::count(cancer_type,mutation_type)
    .data %>>%
      ggplot(aes(x=reorder(mutation_type,mutation_type_order),y=age))+
      geom_violin()+
      geom_boxplot(width=.3,fill="black")+ 
      stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
      facet_wrap(~ cancer_type,ncol=5)+
      ylim(0,120)+
      geom_signif(aes(x=reorder(mutation_type,mutation_type_order),y=age),
                  comparisons = list(c("truncating","none")),
                  test = "wilcox.test",map_signif_level = T, col="red", y_position = 95)+
      geom_signif(aes(x=reorder(mutation_type,mutation_type_order),y=age),
                  comparisons = list(c("damaging","none")),
                  test = "wilcox.test",map_signif_level = T, col="red",y_position = 100)+
      geom_signif(aes(x=reorder(mutation_type,mutation_type_order),y=age),comparisons = list(c("benign","none")),
                  test = "wilcox.test",map_signif_level = T, col="red", y_position = 108)+
      theme(strip.text = element_text(size=6),axis.text.x = element_text(angle = -45, hjust = 0))+
      geom_text(data =.patient_num,aes(x=mutation_type,y=5,label=n),size=3,position="stack")+
      geom_text(data = .mean_age,aes(x=mutation_type,y=120,label=mean),size=3,position = "stack",color="red")
  }

  patient_list %>>%
    filter(!is.na(age)) %>>%
    mutate(patient_focal="focal") %>>%
    left_join(driver_genes %>>%filter(role=="TSG") %>>%dplyr::rename(gene_symbol=gene)%>>%
                mutate(patient_focal="focal") %>>%dplyr::select(gene_symbol,patient_focal))%>>%
    dplyr::select(-patient_focal) %>>%
    left_join(somatic_ranking %>>%filter(mutype=="missense")) %>>%filter(somatic_rank <=10) %>>%
    #left_join(cna_del_rank) %>>%filter(loh_rank <= 30) %>>%
    left_join(nofilter) %>>%
    mutate(mutation_type = ifelse(is.na(mutation_type),"none",mutation_type),
           age=round(age/365.25*100)/100) %>>%
    mutate(mutation_type_order = ifelse(mutation_type == "none", 7, mutation_type_order)) %>>%
    group_by(patient_id)%>>%
    filter(mutation_type_order == min(mutation_type_order)) %>>%
    summarise_all(first) %>>% (?. %>>%count(mutation_type)) %>>%
    class23_ggplot() %>>%
    {ggsave(paste0("age_plot/",.title,".pdf"), width = 15, height = 15)}
}

### all class3 missense or truncate
classed_site %>>% 
  filter(AF_class==3)%>>%
  filter(mutype=="truncating"|mutype=="missense") %>>%
  class23_plot("class3_filter/somatic10")
  class23_plot("class3_missense_or_truncate")

  
### all class2,3 missense_or_truncate
classed_site %>>% 
  filter(AF_class==2|AF_class==3)%>>%
  filter(mutype=="truncating"|mutype=="missense") %>>%
  class23_plot("class23_missense_or_truncate_loh20")
  
  
##################################################################################################
####################################### siteごとに発症年齢分布 ###################################
##################################################################################################

###  class 1 or 2 ###
#all cancer type together
.plot = classed_site %>>%
  filter(mutype!="silent",AF_class == 2,role == "TSG") %>>%
# filter(mutype!="silent",(AF_class == 1 | AF_calss == 2),role == "TSG") %>>%
  mutate(patient_focal="focal") %>>%
  left_join(patient_list %>>%mutate(patient_focal="focal")) %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-AF_class,-patient_focal) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no"))) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt")),
         PolyPhen=ifelse(is.na(PolyPhen),"none",str_extract(PolyPhen,"\\w+")),
         age=round(age/365.25*100)/100) %>>%
  group_by(chr,start,ref,alt,genotype,ac_exac,an_exac) %>>%
  summarise(age = mean(age)) %>>%
  spread(genotype,age) %>>%
  gather(hetero_or_homo,age,refalt,altalt) %>>%
  mutate(diagnosised_age_compared_to_refref = age/refref) %>>%
  ggplot()+
  geom_point(aes(x=ac_exac/an_exac, y=diagnosised_age_compared_to_refref,
                 shape=hetero_or_homo,colour=hetero_or_homo))
.plot
ggsave("age_plot/class2_signtest.pdf",.plot,width = 10,height = 10)

#by cancer type
.plot = classed_site %>>%
  filter(mutype!="silent",AF_class == 2,role == "TSG") %>>%
  mutate(patient_focal="focal") %>>%
  left_join(patient_list %>>%mutate(patient_focal="focal")) %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-AF_class,-patient_focal) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no"))) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt")),
         PolyPhen=ifelse(is.na(PolyPhen),"none",str_extract(PolyPhen,"\\w+")),
         age=round(age/365.25*100)/100) %>>%
  group_by(chr,start,ref,alt,genotype,ac_exac,an_exac,cancer_type) %>>%
  summarise(age = mean(age)) %>>%
  spread(genotype,age) %>>%
  gather(hetero_or_homo,age,refalt,altalt) %>>%
  mutate(diagnosised_age_compared_to_refref = age/refref) %>>%
  ggplot()+
  geom_point(aes(x=ac_exac/an_exac, y=diagnosised_age_compared_to_refref,
                 shape=hetero_or_homo,colour=hetero_or_homo))+
  facet_wrap(~ cancer_type,ncol=5)
.plot
ggsave("age_plot/class2_signtest_byCT.pdf",.plot,height = 15,width = 15)

### class 3 ###
# all cancer type together
.plot = classed_site %>>%
  filter(mutype!="silent",AF_class == 3,role == "TSG") %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-AF_class,-ac_exac,-an_exac) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no")) %>>% dplyr::rename(alt=n_allele2)) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt")),
         age=round(age/365.25*100)/100) %>>%
  group_by(chr,start,ref,alt) %>>%
  summarise(age = mean(age),patient_num=n()) %>>%
  ggplot(aes(x=patient_num, y=age))+
  geom_point()+
  geom_hline(yintercept = round(mean(patient_list$age)/365.25*100)/100)+
  geom_smooth(method = "lm")
.plot
ggsave("age_plot/class3_signtest.pdf",.plot,width = 10,height = 10)

#by cancer type
mean_age = patient_list %>>%group_by(cancer_type) %>>%summarise(age=round(mean(age)/365.25*100)/100)
.plot = classed_site %>>%
  filter(mutype!="silent",AF_class == 3,role == "TSG") %>>%
  dplyr::select(-ac_cancer,-homo_cancer,-AF_class,-ac_exac,-an_exac) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no")) %>>% dplyr::rename(alt=n_allele2)) %>>%
  mutate(genotype = ifelse(is.na(n_genotype),"refref",ifelse(n_genotype=="hetero","refalt","altalt"))) %>>%
  group_by(chr,start,ref,alt) %>>%
  summarise(patient_num=n()) %>>%
  left_join(norm_maf %>>%filter(!(soma_or_germ=="somatic" & LOH =="no")) %>>% dplyr::rename(alt=n_allele2)) %>>%
  group_by(chr,start,ref,alt,patient_num,cancer_type) %>>%
  summarise(age = round(mean(age)/365.25*100)/100) %>>%
  ggplot(aes(x=patient_num, y=age))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~ cancer_type,ncol =5)+
  geom_hline(data = mean_age, aes(yintercept = age))
.plot
ggsave("age_plot/class3_signtest_byCT.pdf",.plot,height = 15,width = 15)


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

##missense の数
missense_count = classed_site %>>%
  filter(AF_class ==3,mutype=="missense") %>>%
  rename(n_allele2=alt) %>>%
  left_join(norm_maf %>>%filter(mutype=="missense",!(soma_or_germ =="somatic" & LOH=="no"),
                                cancer_type!="READ",cancer_type!="KICH")) %>>%
  #left_join(cna_del_rank) %>>%filter(del_rank <40) %>>%
  group_by(cancer_type,patient_id) %>>%
  summarise(missense_num=n()) %>>%
  mutate(missense_count=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
         missense_count_order=ifelse(missense_num >=10,10,missense_num))

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
  ylim(0,90)+
  geom_text(data =missense_count %>>%count(cancer_type,missense_count),
            aes(x=missense_count,y=5,label=n),size=3,position="stack")+
  theme(strip.text = element_text(size=6),axis.text.x = element_text(angle = -45, hjust = 0))
.plot
ggsave("age_plot/cumulative/allgene_missense.pdf",.plot,height = 15,width = 15)


###################################################################################################################
setwd("~/Dropbox/install/tvz/temporary/")
#write_df(somatic_maf,"somatic.maf.gz")
#write_df(varscan_error,"varscan_error")
#write_df(norm_maf,"norm.maf.gz")
#write_df(vcf_exac,"vcf_exac.gz")
#write_df(annotate_ascat,"ascat.gz")
#write_df(driver_genes,"driver_genes.tsv")
#write_df(topdriver_bed,"driver.bed")
#write_df(mid_af_coverage,"AF_mid_coverage_by_patient.tsv.gz")
norm_maf = read_tsv("norm_maf.gz")
varscan_error = read_tsv("varscan_error")
somatic_maf = read_tsv("somatic.maf.gz")
vcf_exac = read_tsv("vcf_exac.gz")
annotate_ascat = read_tsv("ascat.gz")
driver_genes = read_tsv("driver_genes.tsv")
topdriver_bed = read_tsv("driver.bed")
mid_af_coverage = read_tsv("AF_mid_coverage_by_patient.tsv.gz")

