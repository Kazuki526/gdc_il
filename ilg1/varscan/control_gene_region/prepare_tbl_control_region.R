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
control_genes = read_tsv("/Volumes/areca42TB/GRCh38_singlefasta/control_genes.tsv") %>>%
  filter(gene_symbol != "OR8U1") %>>% mutate(focal="yes")
#############################################################################################################
##################################################  nom_maf  ################################################
patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
if(0){
  options(scipen = 10) #startなどを指数表記のまま保存しないため
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
    maf=tibble::tibble(file=list.files(paste0("/Volumes/TLP02-backup2/gdc_download/",.bp,"/maf"))) %>>%
      mutate(filename = paste0('/Volumes/TLP02-backup2/gdc_download/',.bp,"/maf/",file),
             patient_id = str_replace(file,".maf",""))%>>%
      mutate(purrr::map(filename,~strip_maf(.))) %>>%
      unnest() %>>%#(?.)%>>%
      dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                    ref=Reference_Allele,t_allele1=Tumor_Seq_Allele1,t_allele2=Tumor_Seq_Allele2,
                    n_allele1=Match_Norm_Seq_Allele1,n_allele2=Match_Norm_Seq_Allele2) %>>%
      left_join(control_genes) %>>%filter(!is.na(focal)) %>>%
      mutate(soma_or_germ=ifelse((t_allele1==n_allele1)&(t_allele2==n_allele2),"germline","somatic"),
             t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
             n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
      mutate(LOH=ifelse((t_genotype=="homo")&(n_genotype=="hetero"),"LOH","no"))
  }
  tibble::tibble(cancer_type=c("brca","crc","gbm","hnsc","kcc","lgg",
                               "luad","lusc","ov","prad","thca","ucec")) %>>%
    mutate(purrr::map(cancer_type,~extract_norm_maf(.))) %>>%
    unnest() %>>%
    dplyr::select(-file,-filename) %>>%
    mutate(cancer_type = toupper(cancer_type)) %>>%
    write_df("/Volumes/areca42TB2/gdc/control_region/all_patient/control_region.maf.gz")
}

norm_maf_all_cont = read_tsv("/Volumes/areca42TB/tcga/all_patient/control_region.maf.gz") %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  left_join(patient_list) %>>%
  filter(gene_symbol != "GDF2", gene_symbol != "GPRIN2", gene_symbol != "DEFB110") %>>%
  mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH))

######  coverage file
coverage_all = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/coverage_all_cont.tsv.gz")
coverage_male_x = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/coverage_X_male_cont.tsv.gz")

tally_norm_maf_cont = norm_maf_all_cont %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(!(chr=="chrX" & gender=="male" & allele=="n_allele1")) %>>%
  mutate(homo=ifelse(n_genotype=="homo",ifelse(chr=="X" & gender=="male",0,1),0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
            cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
  ungroup() #%>>%
#  left_join(coverage_all_cont) %>>%left_join(coverage_male_x) %>>%
#  mutate(an_male_cancer = ifelse(is.na(an_male_cancer),0,an_male_cancer)) %>>%
#  mutate(an_cancer = an_cancer - an_male_cancer) %>>%dplyr::select(-an_male_cancer)
#############################################################################################################
###############################################  ExAC nonTCGA  ##############################################
strip_maf = function(infile) {
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2',
                'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position",'Transcript_ID')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c', 'c','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
    filter(mutype!="flank",mutype!="splice_region",mutype!="NMD") %>>%
    filter(Consequence!="intron_variant") %>>%
    filter(Consequence!="intron_variant,non_coding_transcript_variant")
}

vcf_exac_cont=read_tsv("/Volumes/areca42TB/exac/control_region/nontcga_liftovered_checked_likevcf.tsv.gz",
                       col_types = "cdccdddddddddddddddddddddddddddddddddd") %>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het))
exac_nonindel_cont=strip_maf("/Volumes/areca42TB/exac/control_region/exac_nontcga_control_region.maf") 
vcf_exac_indel_cont=read_tsv("/Volumes/areca42TB/exac/control_region/nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                             col_types = "cdccdddddddddddddddddddddddddddddddddd")　%>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het)) 
exac_indel_cont = strip_maf("/Volumes/areca42TB/exac/control_region/exac_nontcga_control_region_indel.maf") 

vcf_exac_cont =rbind(vcf_exac_cont,vcf_exac_indel_cont) %>>%
  mutate(chr=paste0("chr",chr))

exac_cont = rbind(exac_nonindel_cont,exac_indel_cont) %>>%left_join(vcf_exac_cont) %>>%
  left_join(control_genes) %>>%filter(!is.na(focal)) %>>%dplyr::select(-focal)
rm(vcf_exac_indel_cont,exac_nonindel_cont,exac_indel_cont)

if(0){
  ###ref_minor_list
  exac_cont %>>%
    filter(AC_Adj/AN_Adj > 0.5| AC_AFR/AN_AFR > 0.5| (AC_FIN+AC_NFE)/(AN_FIN+AN_NFE) > 0.5) %>>%
    write_df("/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_list.tsv")
}
ref_minor_focal_cont = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_coverage_by_patient.tsv.gz")
################################################################################################################
################################################################################################################
################################################ quality filter ################################################
################################################################################################################
################################################################################################################
### HWE test ###
HWE_test_heterom = function(ac,an,hom){
  AF=ac/an
  .matrix = matrix(c(round(an*AF*(1-AF)),round(an/2*AF^2),ac-hom*2,hom), nrow = 2)
  fisher.test(.matrix,alternative = "less")$p.value
}
duplicate_site_cont =exac_cont %>>%
  dplyr::select(gene_symbol,chr,start,ref,alt,mutype,AC_Adj,AN_Adj,AC_Hom) %>>%
  dplyr::rename(ac_exac=AC_Adj,an_exac=AN_Adj,hom_exac=AC_Hom) %>>%
  mutate(hom_exac_hwe=(ac_exac^2/an_exac)/2) %>>%
  filter(chr!="chrX",ac_exac!=0,!is.na(hom_exac)) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE_exac = purrr::map(data,~HWE_test_heterom(.$ac_exac,.$an_exac,.$hom_exac))) %>>%
  unnest() %>>%
  mutate(FDR_exac=p.adjust(HWE_exac)) %>>%filter(FDR_exac<0.01) %>>%
#  full_join(tally_norm_maf%>>%
#              dplyr::select(chr,start,ref,alt,ac_cancer,hom_cancer,gene_symbol,mutype) %>>%
#              left_join(coverage_all ) %>>%
#              mutate(hom_cancer_hwe = ac_cancer^2/an_cancer/2) %>>%
#              filter(chr!="chrX") %>>%
#              nest(-chr,-start,-ref,-alt) %>>%
#              mutate(HWE_cancer=purrr::map(data,~HWE_test_heterom(.$ac_cancer,.$an_cancer,.$hom_cancer))) %>>%
#              unnest() %>>%
#              mutate(FDR_cancer=p.adjust(HWE_cancer)) %>>%filter(FDR_cancer<0.01)) %>>%
  filter(ref!="-" & alt!="-")
##### somaticでrecurrentなmutationはgermで起こっているとは考えにくい（様々なエラーが考えられる）
#いちおうEXACでAF>1%となっているsiteはこれに含めない
somatic_recurrent_cont = norm_maf_all_cont%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)

#normalでaltalt, tumorでrefaltとなってる際にnormalでrefのdepth=0のものだけ採用！
#また同じサイトでこのエラーが有意に多い(100patient以上の)siteは解析に使用しないことにした。(3site)
#しかし control_regionでは存在しない、、
varscan_error_cont = norm_maf_all_cont %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(an_error = n()*2) %>>%
  ungroup() 
varscan_error_site_cont = norm_maf_all_cont %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  dplyr::select(patient_id,chr,start,ref,alt) %>>%
  mutate(varscan_error_focal="yes")


quality_filter_cont =
  function(.data,.data_type="vcf",.fdr=0.01,.database="all",.duplicate=T,.somatic=T,.varscan=F){
  .site = duplicate_site_cont
  if(.database!="all"){
    if(.database=="cancer"){
      .site = .site %>>% filter(FDR_cancer < .fdr)
    }else if(.database=="exac"){
      .site = .site %>>% filter(FDR_exac < .fdr)
    }else{stop(paste0("database variable is wrong .database=",.database,
                      "\ncancer, exac is correct."))}
  }
  if(.varscan){
    if(.data_type=="vcf"){
      .data = .data %>>%
        left_join(varscan_error_cont) %>>%
        mutate(ac_cancer = ifelse(!is.na(an_cancer),ac_cancer-an_error,ac_cancer),
               an_cancer = ifelse(!is.na(an_cancer),an_cancer-an_error,an_cancer) - an_male_cancer) %>>%
        dplyr::select(-an_male_cancer)
    }else if(.data_type=="maf"){
      .data = .data %>>%
        left_join(varscan_error_site_cont) %>>%
        filter(is.na(varscan_error_focal)) %>>%
        group_by(patient_id,gene_symbol,Protein_position) %>>%
        mutate(same_codon=n()) %>>%
        filter(same_codon==1 & !str_detect(Protein_position,"-")) %>>%
        ungroup() %>>%
        dplyr::select(-varscan_error_focal, -same_codon)
    }else{stop(paste0(".datatype variabel is wrong .datatype=",.type,"\nvcf, maf is correct"))}
  }
  if(.data_type=="maf"){
    .data = .data %>>%
      filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
      mutate(alt=n_allele2)
  }
  kb_sagittal = function(.start){
    .head=.start-1001
    .tail=.start+1000
    data.frame(start=.head:.tail)
  }
  .remove_site = .site %>>%
    dplyr::select(chr,start) %>>%
    mutate(start = purrr::map(start,~kb_sagittal(.)),duplicate_focal = "yes") %>>%
    unnest() %>>%distinct()
  .data = .data %>>%
    left_join(.remove_site) %>>%
    left_join(somatic_recurrent_cont) %>>%
    filter(if(.duplicate==T){is.na(duplicate_focal)}else{chr==chr}) %>>%
    filter(if(.somatic==T){is.na(recurrent_focal)}else{chr==chr})%>>%
    dplyr::select(-recurrent_focal,-duplicate_focal)
}


