######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.Rを実行してからこのスクリプトを開始する！！！

#prepare 1000genomes UK10K EXAC(nonTCGA) data sets
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
###1kg###

sample_num_1kg=5008
sample_num_X_1kg=sample_num_1kg - 1233
#vcf_1kg_ac0=read_tsv("/working/1000genomes/maf/extract/ac0_variation.tsv")
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",col_types = "ccccddcccccccccc") #germ_pvalue_plot.Rにて作成

super_population = function(.data){
  mutate(.data, SPC = dplyr::recode(Population, #SPC = Super Population Code
                                    CHB = "EAS", JPT = "EAS", CHS = "EAS", CDX = "EAS", KHV = "EAS",
                                    CEU = "EUR", TSI = "EUR", FIN = "EUR", GBR = "EUR", IBS= "EUR",
                                    YRI = "AFR", LWK = "AFR", GWD = "AFR", MSL = "AFR", ESN = "AFR",
                                    ASW = "AFR", ACB = "AFR", MXL = "AMR", PUR = "AMR", CLM = "AMR",
                                    PEL = "AMR", GIH = "SAS", PJL = "SAS", BEB = "SAS", STU = "SAS",
                                    ITU = "SAS"))
} ## AFR = African, AMR = Ad Mixed American, EAS = East Asian, EUR = European, SAS = South Asian
thousand_genome_info=read_csv("~/git/all_patient/1kg_sample_info.csv") %>>%
  super_population() %>>%dplyr::rename(sample = Sample) %>>%
  dplyr::select(sample,SPC)

tally_1kg_all= maf_1kg_all %>>%
  left_join(thousand_genome_info) %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_1kg=n()) %>>% ungroup() %>>%
  left_join(maf_1kg_all %>>% group_by(sample,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(hom_1kg=n()/2)) %>>%
  mutate(an_1kg = ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg))%>>%
  mutate(chr=paste0("chr",chr),hom_1kg=ifelse(is.na(hom_1kg),0,hom_1kg)) #%>>%full_join(vcf_1kg_ac0)

rm(sample_num_1kg,sample_num_X_1kg,maf_1kg_all)
####UK 10K ###
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle_likevcf.tsv.gz",col_types = "cdccdddd")
uk10k=strip_maf("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle.maf") %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand))%>>%
  left_join(vcf_10k%>>%mutate(start=ifelse(ref=="-",start -1,start))) %>>%
  mutate(chr=paste0("chr",chr)) %>>%
  dplyr::rename(ac_uk=uk_ac,an_uk=uk_an,hom_uk=uk_althomo,het_uk=uk_hetero) %>>%
  filter(an_uk >3000)
rm(vcf_10k)

### EXAC nonTCGA ###
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/nontcga_liftovered_checked_likevcf.tsv.gz",
                  col_types = "cdccdddddddddddddddddddddddddddddddddd") %>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het))
exac_nonindel=strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver.maf") %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand))
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                        col_types = "cdccdddddddddddddddddddddddddddddddddd")　%>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het)) 
exac_indel = strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver_indel.maf") %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand))

vcf_exac =rbind(vcf_exac,vcf_exac_indel) %>>%
  mutate(chr=paste0("chr",chr))

exac = rbind(exac_nonindel,exac_indel) %>>%left_join(vcf_exac)
rm(vcf_exac_indel,exac_nonindel,exac_indel)


#######################################################################################################
################################filtering duplicate? position #########################################
#######################################################################################################
### HWE test ###
if(0){
HWE_test_heterom = function(ac,an,hom){
  AF=ac/an
  .matrix = matrix(c(round(an*AF*(1-AF)),round(an/2*AF^2),ac-hom*2,hom), nrow = 2)
  fisher.test(.matrix,alternative = "less")$p.value
}

duplicate_site = exac %>>%
  dplyr::select(gene_symbol,chr,start,ref,alt,mutype,AC_Adj,AN_Adj,AC_Hom) %>>%
  dplyr::rename(ac_exac=AC_Adj,an_exac=AN_Adj,hom_exac=AC_Hom) %>>%
  mutate(hom_exac_hwe=(ac_exac^2/an_exac)/2) %>>%
  filter(chr!="chrX",ac_exac!=0,!is.na(hom_exac)) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE_exac = purrr::map(data,~HWE_test_heterom(.$ac_exac,.$an_exac,.$hom_exac))) %>>%
  unnest() %>>%
  mutate(FDR_exac=p.adjust(HWE_exac)) %>>%filter(FDR_exac<0.01)%>>%
  full_join(tally_norm_maf%>>%
              dplyr::select(chr,start,ref,alt,ac_cancer,hom_cancer,gene_symbol,mutype,an_cancer) %>>%
              mutate(hom_cancer_hwe = ac_cancer^2/an_cancer/2) %>>%
              filter(chr!="chrX") %>>%
              nest(-chr,-start,-ref,-alt) %>>%
              mutate(HWE_cancer=purrr::map(data,~HWE_test_heterom(.$ac_cancer,.$an_cancer,.$hom_cancer))) %>>%
              unnest() %>>%
              mutate(FDR_cancer=p.adjust(HWE_cancer)) %>>%filter(FDR_cancer<0.01)) %>>%
  full_join(tally_1kg_all %>>%dplyr::select(chr,start,ref,alt,ac_1kg,an_1kg,hom_1kg,gene_symbol,mutype) %>>%
              mutate(hom_1kg_hwe = ac_1kg^2/an_1kg/2) %>>%
              filter(chr!="chrX") %>>%
              nest(-chr,-start,-ref,-alt) %>>%
              mutate(HWE_1kg = purrr::map(data,~HWE_test_heterom(.$ac_1kg,.$an_1kg,.$hom_1kg)))%>>%
              unnest() %>>%
              mutate(FDR_1kg=p.adjust(HWE_1kg)) %>>%filter(FDR_1kg<0.01)) %>>%
  full_join(uk10k %>>%dplyr::select(chr,start,ref,alt,ac_uk,an_uk,hom_uk,gene_symbol,mutype) %>>%
              mutate(hom_uk_hwe = ac_uk^2/an_uk/2) %>>%
              filter(chr!="chrX") %>>%
              nest(-chr,-start,-ref,-alt) %>>%
              mutate(HWE_uk = purrr::map(data,~HWE_test_heterom(.$ac_uk,.$an_uk,.$hom_uk))) %>>%
              unnest() %>>%
              mutate(FDR_uk=p.adjust(HWE_uk)) %>>%filter(FDR_uk<0.01))%>>%
  filter(ref!="-" & alt!="-")
write_df(duplicate_site,"/Volumes/areca42TB2/gdc/varscan/all_patient/duplicate_site.tsv")
}
duplicate_site = read_tsv("/Volumes/areca42TB2/gdc/varscan/all_patient/duplicate_site.tsv")

##### somaticでrecurrentなmutationはgermで起こっているとは考えにくい（様々なエラーが考えられる）
#いちおうEXACでAF>1%となっているsiteはこれに含めない
somatic_recurrent = norm_maf_all%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)

#normalでaltalt, tumorでrefaltとなってる際にnormalでrefのdepth=0のものだけ採用！
#また同じサイトでこのエラーが有意に多い(100patient以上の)siteは解析に使用しないことにした。(3site)
varscan_error = norm_maf_all %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(an_error = n()*2) %>>%
  ungroup() 
varscan_error_site = norm_maf_all %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  dplyr::select(patient_id,chr,start,ref,alt) %>>%
  mutate(varscan_error_focal="yes")


quality_filter= function(.data,.data_type="vcf",.fdr=0.01,.database="all",.duplicate=T,.somatic=T,.varscan=F){
  .site = duplicate_site
  if(.database!="all"){
    if(.database=="cancer"){
      .site = .site %>>% filter(FDR_cancer < .fdr)
    }else if(.database=="exac"){
      .site = .site %>>% filter(FDR_exac < .fdr)
    }else if(.database=="1kg"){
      .site = .site %>>% filter(FDR_1kg < .fdr)
    }else if(.database=="uk"){
      .site = .site %>>% filter(FDR_uk < .fdr)
    }else{stop(paste0("database variable is wrong .database=",.database,
                      "\ncancer, exac, 1kg, uk is correct."))}
  }
  if(.varscan){
    if(.data_type=="vcf"){
      .data = .data %>>%
        left_join(varscan_error) %>>%
        mutate(ac_cancer = ifelse(!is.na(an_cancer),ac_cancer-an_error,ac_cancer),
               an_cancer = ifelse(!is.na(an_cancer),an_cancer-an_error,an_cancer) - an_male_cancer) %>>%
        dplyr::select(-an_male_cancer)
    }else if(.data_type=="maf"){
      .data = .data %>>%
        left_join(varscan_error_site) %>>%
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
    filter(gene_symbol!="KMT2C") %>>%
    dplyr::select(chr,start) %>>%
    mutate(start = purrr::map(start,~kb_sagittal(.)),duplicate_focal = "yes") %>>%
    unnest() %>>%distinct()
  .data = .data %>>%
    left_join(.remove_site) %>>%
    left_join(somatic_recurrent) %>>%
    filter(if(.duplicate){gene_symbol!="KMT2C"}else{chr==chr}) %>>%
    filter(if(.duplicate==T){is.na(duplicate_focal)}else{chr==chr}) %>>%
    filter(if(.somatic==T){is.na(recurrent_focal)}else{chr==chr})%>>%
    dplyr::select(-recurrent_focal,-duplicate_focal)
}
s                                                          
