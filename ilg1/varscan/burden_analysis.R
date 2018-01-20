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
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand))%>>%
    filter(mutype!="flank",mutype!="splice_region") %>>%
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
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_1kg=n())) %>>%
  mutate(an_1kg = ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg))%>>%
  mutate(chr=paste0("chr",chr),homo_1kg=ifelse(is.na(homo_1kg),0,homo_1kg)) #%>>%
  #full_join(vcf_1kg_ac0)

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
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",
                  col_types = "cdccdddddddddddddddddddddddddddddddddd") 
exac_nonindel=strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver.maf") 
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                        col_types = "cdccdddddddddddddddddddddddddddddddddd")
exac_indel = strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver_indel.maf") 

vcf_exac =rbind(vcf_exac,vcf_exac_indel) %>>%
  mutate(chr=paste0("chr",chr))
  
exac = rbind(exac_nonindel,exac_indel) %>>%left_join(vcf_exac)
rm(vcf_exac_indel,exac_nonindel,exac_indel)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#ここから解析
#kokokarakaiseki
#here start analysis
#####################################################################################################################
### Charlesら(NatCom2015の解析の確認)
####### white ######
##################################
### HDD返ってくるまでの緊急用
white_prop=4424/6514
black_prop=1412/6514
cancer_type_prop = norm_maf_all %>>%
  count(cancer_type,patient_id,age)%>>%
  left_join(patient_race)%>>%count(cancer_type,race) %>>%
  mutate(prop=nn/6514)
##################################
# splice & truncating (ref minor でtruncating or spliceはないのでとりあえず無視)
truncating_focal_site = quality_filter(norm_maf_all,.data_type="maf") %>>%
  filter(mutype == "splice"| mutype =="truncating") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  left_join(coverage_all) %>>%
  full_join(quality_filter(exac)) %>>% filter(mutype == "splice"| mutype =="truncating") %>>%
  filter((is.na(AC_Adj)|AC_Adj/AN_Adj *100 <0.05), sum(AC_Adj,ac_cancer,na.rm = T) >2, chr!="chrX") %>>%
  mutate(focal = "ok") %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,focal)

cancer_white = quality_filter(norm_maf_all,.data_type="maf") %>>%
  left_join(patient_race) %>>%filter(race=="white") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt, chr !="chrX")%>>%
  group_by(cancer_type,chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,cancer_type,ac_cancer)%>>%
  left_join(coverage_all)
exac_white = exac %>>% quality_filter() %>>%
  filter(chr != "chrX") %>>%
  mutate(ac_cont = AC_FIN + AC_NFE,an_cont = AN_FIN + AN_NFE) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont)

#### all cancer type ####
TFT = function(.tbl){
  case_f =.tbl$taf_cancer
  control_f = .tbl$taf_cont
  case =.tbl$an_cancer
  control = .tbl$an_cont
  ifelse(case_f ==0 & control_f==0, 1,
    prop.test(c(case_f*case,control_f*control),c(case,control),alternative = "greater")$p.value)
}
TFT_all_cancer_type = 
  function(.tbl,.cont_tbl,.focal_site=truncating_focal_site,.race=white_prop){
    full_join(.tbl %>>% group_by(chr,start,end,ref,alt) %>>%
                summarise(ac_cancer=sum(ac_cancer),gene_symbol=first(gene_symbol),
                          Consequence=first(Consequence),mutype=first(mutype),an_cancer=first(an_cancer)),
              .cont_tbl) %>>%
      left_join(.focal_site) %>>%filter(!is.na(focal)) %>>%
      mutate(an_cancer = ifelse(is.na(an_cancer),1,round(an_cancer*.race)), ####################################
             an_cont   = ifelse(is.na(an_cont  ),1,an_cont  )) %>>%
      mutate(af_cancer = ifelse(is.na(ac_cancer),0,ac_cancer/an_cancer),
             af_cont   = ifelse(is.na(ac_cont  ),0,ac_cont  /an_cont  )) %>>%
      group_by(gene_symbol) %>>%
      summarise(taf_cancer = sum(af_cancer,na.rm = T),taf_cont =sum(af_cont,na.rm = T),
                an_cancer = max(an_cancer,na.rm = T),an_cont = max(an_cont,na.rm = T)) %>>%
      nest(-gene_symbol)%>>%
      mutate(tft = purrr::map(data, ~TFT(.))) %>>%unnest()%>>%
      {left_join(driver_genes %>>%dplyr::rename(gene_symbol=gene)%>>%dplyr::select(gene_symbol,role),.)} %>>%
      mutate(tft =ifelse(is.na(tft),1,tft)) %>>%
      mutate(FDR = p.adjust(tft,"fdr"))
  }
total_frequency_truncate = TFT_all_cancer_type(cancer_white,exac_white)
#### plot ####
FDR_barplot = function(.tbl,.FDR=0.5){
  log_reverse = function(.log){
    parse(text = paste0("10^",-log10(.log)))
  }
  .tbl %>>% filter(FDR < .FDR)%>>%
    arrange(FDR) %>>%
    mutate(FDR=FDR*(100^(-log10(FDR))))%>>%
    ggplot(aes(x=reorder(gene_symbol,desc(FDR)),y=FDR,fill=gene_symbol))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 50,size=0.2)+
    geom_hline(yintercept = 0,size=1,colour="gray")+
    scale_y_log10(labels = log_reverse)+
    #scale_fill_brewer(palette="Set1")+
    guides(fill="none")+
    theme_bw()+
    theme( panel.grid.major.x = element_blank(),
           panel.background = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
}
.plot=FDR_barplot(total_frequency_truncate)
.plot
ggsave("burden_plot/truncate_fdr_white_likeCharles.pdf",.plot,height = 7,width =4)
#ggsave("burden_plot/truncate_fdr_likeCharles_cut_singleduble.pdf",.plot,height = 4,width = 2.5)

#### by cancer type ####
TFT_by_cancer_type = function(.tbl,.control,.focal_site=truncating_focal_site){
  full_join(.tbl,.control) %>>%
    left_join(.focal_site) %>>%filter(!is.na(focal)) %>>%
    #filter(mutype == "truncating" | mutype == "splice") %>>% #singleton, doubleton cutしない時用
    mutate(an_cancer = ifelse(is.na(an_cancer),1,an_cancer),
           an_cont   = ifelse(is.na(an_cont  ),1,an_cont  )) %>>%
    mutate(af_cancer = ifelse(is.na(ac_cancer),0,ac_cancer/an_cancer),
           af_cont   = ifelse(is.na(ac_cont  ),0,ac_cont  /an_cont  )) %>>%
    group_by(gene_symbol) %>>%
    summarise(taf_cancer = sum(af_cancer,na.rm = T),taf_cont=sum(af_cont,na.rm = T),
              an_cancer = max(an_cancer,na.rm = T),an_cont = max(an_cont,na.rm = T)) %>>%
    nest(-gene_symbol) %>>%
    mutate(tft = purrr::map(data, ~TFT(.))) %>>%unnest()%>>%
    {left_join(driver_genes %>>%dplyr::rename(gene_symbol=gene)%>>%dplyr::select(gene_symbol,role),.)} %>>%
    mutate(tft =ifelse(is.na(tft),1,tft)) %>>%
    mutate(FDR = p.adjust(tft,"fdr"))
}
total_frequency_truncate_by_cancertype = cancer_white %>>%
  left_join(cancer_type_prop %>>% filter(race == "white")) %>>%
  mutate(an_cancer = round(an_cancer * prop)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type(.,exac_white))) %>>%unnest()

######## bubble plot #########
bubble_plot = function(.tbl_by_CT,.tbl_all_CT,.FDR=0.1){
  rbind(.tbl_by_CT,
        .tbl_all_CT %>>% mutate(cancer_type = "Pan-cancer")) %>>%
    group_by(gene_symbol) %>>%
    mutate(focal = ifelse(all(FDR > .FDR),"no","ok")) %>>% #singleton, doubleton cutしなかった時用
    filter(focal=="ok") %>>%ungroup() %>>%
    mutate(canty_order=ifelse(cancer_type=="Pan-cancer",2,1))%>>%
    ggplot()+
    geom_point(aes(x=reorder(as.factor(cancer_type),canty_order),
                   y=reorder(as.factor(gene_symbol),desc(as.factor(gene_symbol))),
                   size=taf_cancer*100,color=cancer_type))+
    geom_point(data = .tbl_by_CT%>>%filter(FDR < 0.05),
               aes(x=cancer_type,y=gene_symbol,size=taf_cancer*300),shape =21,stroke = 1)+
    #scale_color_brewer(palette = "Set1")+
    guides(colour="none",size= guide_legend("% of samples "))+
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
}
.plot = bubble_plot(total_frequency_truncate_by_cancertype,
              total_frequency_truncate)
.plot
ggsave("burden_plot/truncate_bubble_plot_white.pdf",.plot,height = 7,width = 4)
#ggsave("burden_plot/truncate_bubble_plot_cut_singledouble.pdf",.plot,height = 7,width = 4)
###################################################################################################################
###################################################################################################################
###################################################################################################################
#ここからはサプリ用
# black
cancer_black = quality_filter(norm_maf_all,.data_type="maf") %>>%
  left_join(patient_race) %>>%filter(race=="black") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt, chr !="chrX")%>>%
  group_by(cancer_type,chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,cancer_type,ac_cancer)%>>%
  left_join(coverage_all)
exac_black = exac %>>% quality_filter() %>>%
  filter(chr != "chrX") %>>%
  mutate(ac_cont = AC_AFR, an_cont = AN_AFR) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont)
total_frequency_truncate_black = TFT_all_cancer_type(cancer_black,exac_black,.race = black_prop)
.plot = FDR_barplot(total_frequency_truncate_black,.FDR = 0.5)
.plot
ggsave("burden_plot/truncate_fdr_black_likeCharles.pdf",.plot,height = 7, width = 4)
############FDR有意なものがない！！！？？？###################

total_frequency_truncate_by_cancertype_balck = cancer_black %>>%
  left_join(cancer_type_prop %>>% filter(race == "black")) %>>%
  mutate(an_cancer = round(an_cancer * prop)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type(.,exac_black))) %>>%unnest()
.plot = bubble_plot(total_frequency_truncate_by_cancertype_balck,
                    total_frequency_truncate_black)
.plot
ggsave("burden_plot/truncate_bubble_plot_black.pdf",.plot,height = 7,width = 4)
####################################################################################################################
## 1000 genomes
sample_num_1kg_white=1338
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",col_types = "ccccddcccccccccc") #germ_pvalue_plot.Rにて作成
tally_1kg_white= maf_1kg_all %>>%
  left_join(thousand_genome_info) %>>%filter(SPC=="EUR") %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_cont=n()) %>>%ungroup() %>>%
  mutate(chr=paste0("chr",chr),an_cont=sample_num_1kg_white) %>>%
  quality_filter()
rm(maf_1kg_all)
total_frequency_truncate_1kg = TFT_all_cancer_type(cancer_white,tally_1kg_white)
.plot = FDR_barplot(total_frequency_truncate_1kg,.FDR = 0.6)
.plot
#####またまた 有意なものがない、、、
ggsave("burden_plot/1kg/truncate_fdr_white_likeCharles.pdf",.plot,height = 3, width = 4)

total_frequency_truncate_by_cancertype_1kg = cancer_white %>>%
  left_join(cancer_type_prop %>>% filter(race == "white")) %>>%
  mutate(an_cancer = round(an_cancer * prop)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type(.,tally_1kg_white))) %>>%unnest()
.plot = bubble_plot(total_frequency_truncate_by_cancertype_1kg,
                    total_frequency_truncate_1kg)
.plot
######OVのBRCA1だけ有意
ggsave("burden_plot/1kg/truncate_bubble_plot_white.pdf",.plot,height = 3,width = 4)

####################################################################################################################
## UK10K
truncating_focal_site_uk = quality_filter(norm_maf_all,.data_type="maf") %>>%
  filter(mutype == "splice"| mutype =="truncating") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  full_join(quality_filter(exac) %>>% filter(mutype == "splice"| mutype =="truncating"))%>>%
  full_join(quality_filter(uk10k) %>>%filter(mutype == "splice"| mutype =="truncating"))%>>%
  filter((is.na(AC_Adj)|AC_Adj/AN_Adj *100 <0.05), sum(AC_Adj,ac_cancer,ac_uk,na.rm = T) >2, chr!="chrX") %>>%
  mutate(focal = "ok") %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,focal)
uk_tft = uk10k %>>% quality_filter() %>>%
  dplyr::rename(ac_cont = ac_uk,an_cont = an_uk) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont)
total_frequency_truncate_uk =
  TFT_all_cancer_type(cancer_white, uk_tft,
                      .focal_site = truncating_focal_site_uk)
.plot = FDR_barplot(total_frequency_truncate_uk,.FDR = 0.6)
.plot
ggsave("burden_plot/uk10k/truncate_fdr_white_likeCharles.pdf",.plot,height = 3, width = 4)
total_frequency_truncate_by_cancertype_uk = cancer_white %>>%
  left_join(cancer_type_prop %>>% filter(race == "white")) %>>%
  mutate(an_cancer = round(an_cancer * prop)) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type(.,uk_tft,truncating_focal_site_uk))) %>>%
  unnest()
.plot = bubble_plot(total_frequency_truncate_by_cancertype_uk,
                    total_frequency_truncate_uk)
.plot
ggsave("burden_plot/uk10k/truncate_bubble_plot_white.pdf",.plot,height = 3, width = 4)

