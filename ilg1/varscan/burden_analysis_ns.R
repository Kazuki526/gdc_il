######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.Rを実行してからこのスクリプトを開始する！！！
loadNamespace('cowplot')

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

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### all cancer type ####
TFT = function(.tbl){
  case_f =.tbl$taf_cancer
  control_f = .tbl$taf_cont
  case =.tbl$an_cancer
  control = .tbl$an_cont
  ifelse(case_f ==0 & control_f==0, 1,
         prop.test(c(case_f*case,control_f*control),c(case,control),alternative = "greater")$p.value)
}
TFT_all_cancer_type_ns = 
  function(.tbl,.cont_tbl,.focal_site=ns_focal_site){
    full_join(.tbl %>>% group_by(chr,start,end,ref,alt) %>>%
                summarise(ac_cancer=sum(ac_cancer), an_cancer=sum(an_cancer),
                          an_all_cantype = first(an_all_cantype), gene_symbol=first(gene_symbol),
                          Consequence=first(Consequence), mutype=first(mutype)),
              .cont_tbl) %>>%
      left_join(.focal_site) %>>%filter(!is.na(focal)) %>>%
      mutate(an_all_cantype = ifelse(is.na(an_all_cantype),1,an_all_cantype),
             an_cont   = ifelse(is.na(an_cont  ),1,an_cont  )) %>>%
      mutate(af_cancer = ifelse(is.na(ac_cancer),0,ac_cancer/an_all_cantype),
             af_cont   = ifelse(is.na(ac_cont  ),0,ac_cont  /an_cont  )) %>>%View
      group_by(gene_symbol,mutype) %>>%
      summarise(taf_cancer = sum(af_cancer,na.rm = T),taf_cont =sum(af_cont,na.rm = T),
                an_cancer = max(an_all_cantype,na.rm = T),an_cont = max(an_cont,na.rm = T)) %>>%
      ungroup() %>>%
      nest(-gene_symbol,-mutype)%>>%
      mutate(tft = purrr::map(data, ~TFT(.))) %>>%unnest()%>>%
      {left_join(driver_genes %>>%dplyr::rename(gene_symbol=gene)%>>%dplyr::select(gene_symbol,role),.)} %>>%
      mutate(tft =ifelse(is.na(tft),1,tft)) %>>%
      mutate(FDR = p.adjust(tft,"fdr"))
  }
#### plot ####
FDR_barplot_ns = function(.tbl,.FDR=0.2){
  log_reverse = function(.log){
    parse(text = paste0("10^",-log10(.log)))
  }
  .tbl =.tbl %>>% filter(FDR < .FDR)%>>%
    arrange(FDR) %>>%
    mutate(FDR=10^(-log10(FDR)))
  .plot_fdr = function(.tibble){
    .tibble %>>%
      ggplot(aes(x=reorder(gene_symbol,desc(FDR)),y=FDR,fill=gene_symbol))+
      geom_bar(stat = "identity")+
      geom_hline(yintercept = 20,size=0.2)+
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
  .mis_plot = .plot_fdr(.tbl %>>%filter(mutype=="missense"))+ggtitle("nonsynonymous")
  .sil_plot = .plot_fdr(.tbl %>>%filter(mutype=="silent"))+ggtitle("synonymous")
  .plot = cowplot::plot_grid(.mis_plot,.sil_plot, ncol = 1)
  plot(.plot)
  return(.plot)
}
#### by cancer type ####
TFT_by_cancer_type_ns = function(.tbl,.control,.focal_site=ns_focal_site){
  full_join(.tbl,.control) %>>%
    left_join(.focal_site) %>>%filter(!is.na(focal)) %>>%
    mutate(an_cancer = ifelse(is.na(an_cancer),1,an_cancer),
           an_cont   = ifelse(is.na(an_cont  ),1,an_cont  )) %>>%
    mutate(af_cancer = ifelse(is.na(ac_cancer),0,ac_cancer/an_cancer),
           af_cont   = ifelse(is.na(ac_cont  ),0,ac_cont  /an_cont  )) %>>%
    group_by(gene_symbol,mutype) %>>%
    summarise(taf_cancer = sum(af_cancer,na.rm = T),taf_cont=sum(af_cont,na.rm = T),
              an_cancer = max(an_cancer,na.rm = T),an_cont = max(an_cont,na.rm = T)) %>>%
    ungroup() %>>%
    nest(-gene_symbol,-mutype) %>>%
    mutate(tft = purrr::map(data, ~TFT(.))) %>>%unnest()%>>%
    {left_join(driver_genes %>>%dplyr::rename(gene_symbol=gene)%>>%dplyr::select(gene_symbol,role),.)} %>>%
    mutate(tft =ifelse(is.na(tft),1,tft)) %>>%
    mutate(FDR = p.adjust(tft,"fdr"))
}
######## bubble plot #########
bubble_plot_ns = function(.tbl_by_CT,.tbl_all_CT,.FDR=0.1){
  .tbl = rbind(.tbl_by_CT,
               .tbl_all_CT %>>% mutate(cancer_type = "Pan-cancer")) %>>%
    group_by(gene_symbol,mutype) %>>%
    mutate(focal = ifelse(all(FDR > .FDR),"no","ok")) %>>% #singleton, doubleton cutしなかった時用
    filter(focal=="ok") %>>%ungroup() %>>%
    mutate(canty_order=ifelse(cancer_type=="Pan-cancer",2,1)) %>>%
    filter(taf_cancer != 0)
  .plot_bubble = function(.tibble){
    .tibble %>>%
      ggplot()+
      geom_point(aes(x=reorder(as.factor(cancer_type),canty_order),
                     y=reorder(as.factor(gene_symbol),desc(as.factor(gene_symbol))),
                     size=taf_cancer*100,color=cancer_type))+
      geom_point(data = .tibble%>>%filter(FDR < 0.05),
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
  .mis_bubble = .plot_bubble(.tbl %>>% filter(mutype == "missense"))+ggtitle("nonsynonymous")
  .sil_bubble = .plot_bubble(.tbl %>>% filter(mutype == "silent"))+ggtitle("synonymous")
  .plot = cowplot::plot_grid(.mis_bubble,.sil_bubble, ncol = 2)
  plot(.plot)
  return(.plot)
}
##################################
# white
ns_focal_site = quality_filter(norm_maf_all,.data_type="maf") %>>%
  filter(mutype == "missense"| mutype =="silent") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  left_join(coverage_all) %>>%
  full_join(quality_filter(exac)) %>>% filter(mutype == "missense"| mutype =="silent") %>>%
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
  left_join(coverage_all_by_cancer_type %>>%
              dplyr::rename(an_cancer =an_white) %>>%
              dplyr::select(-an_black,-an_other)) %>>%
  left_join(coverage_all_by_cancer_type %>>%
              group_by(chr,start) %>>%summarise(an_all_cantype = sum(an_white)) %>>%
              ungroup()) %>>%
  filter(an_all_cantype > max(an_all_cantype)*0.5)
exac_white = exac %>>% quality_filter() %>>%
  filter(chr != "chrX") %>>%
  mutate(ac_cont = AC_FIN + AC_NFE,an_cont = AN_FIN + AN_NFE) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont)%>>%
  filter(an_cont > max(an_cont)*0.5)

total_frequency_ns = TFT_all_cancer_type_ns(cancer_white,exac_white)
.plot_tft_w = FDR_barplot_ns(total_frequency_ns,.FDR = 0.1)
ggsave("burden_plot/ns_fdr_white_likeCharles.pdf",.plot_tft_w,height = 10,width =4)
#ggsave("burden_plot/truncate_fdr_likeCharles_cut_singleduble.pdf",.plot,height = 4,width = 2.5)

total_frequency_ns_by_cancertype = cancer_white %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type_ns(.,exac_white))) %>>%unnest()
.plot_bb_w = bubble_plot_ns(total_frequency_ns_by_cancertype,total_frequency_ns,.FDR = 0.05)
ggsave("burden_plot/ns_bubble_plot_white.pdf",.plot_bb_w,height = 7,width = 8)
############################################
#black
ns_focal_site_black = quality_filter(norm_maf_all,.data_type="maf") %>>%
  filter(mutype == "missense"| mutype =="silent") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  left_join(coverage_all) %>>%
  full_join(quality_filter(exac)) %>>% filter(mutype == "missense"| mutype =="silent") %>>%
  filter((is.na(AC_AFR)|AC_AFR/AN_AFR *100 <0.05), sum(AC_Adj,ac_cancer,na.rm = T) >2, chr!="chrX") %>>%
  mutate(focal = "ok") %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,focal) 
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
  left_join(coverage_all_by_cancer_type %>>%
              dplyr::rename(an_cancer =an_black) %>>%
              dplyr::select(-an_white,-an_other)) %>>%
  left_join(coverage_all_by_cancer_type %>>%
              group_by(chr,start) %>>%summarise(an_all_cantype = sum(an_black)) %>>%
              ungroup())%>>%
  filter(an_all_cantype > max(an_all_cantype)*0.5)
exac_black = exac %>>% quality_filter() %>>%
  filter(chr != "chrX") %>>%
  mutate(ac_cont = AC_AFR, an_cont = AN_AFR) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont)%>>%
  filter(an_cont > max(an_cont)*0.5)
total_frequency_ns_black = TFT_all_cancer_type_ns(cancer_black,exac_black,
                                                  .focal_site = ns_focal_site_black)
.plot_tft_b = FDR_barplot_ns(total_frequency_ns_black,.FDR = 0.2)
ggsave("burden_plot/ns_fdr_black_likeCharles.pdf",.plot_tft_b,height = 7, width = 2)


total_frequency_ns_by_cancertype_black = cancer_black %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type_ns(
    .,exac_black,.focal_site=ns_focal_site_black))) %>>%unnest()
.plot_bb_b = bubble_plot_ns(total_frequency_ns_by_cancertype_black,total_frequency_ns_black)
ggsave("burden_plot/ns_bubble_plot_black.pdf",.plot_bb_b,height = 5,width = 8)

## 1000 genomes
sample_num_1kg_white=1338
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",
                     col_types = "ccccddcccccccccc") #germ_pvalue_plot.Rにて作成
tally_1kg_white= maf_1kg_all %>>%
  left_join(thousand_genome_info) %>>%filter(SPC=="EUR") %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_cont=n()) %>>%ungroup() %>>%
  mutate(chr=paste0("chr",chr),an_cont=sample_num_1kg_white) %>>%
  quality_filter()
rm(maf_1kg_all)
total_frequency_ns_1kg = TFT_all_cancer_type_ns(cancer_white,tally_1kg_white)
.plot_tft_1 = FDR_barplot_ns(total_frequency_ns_1kg,.FDR = 0.2)
ggsave("burden_plot/1kg/ns_fdr_white_likeCharles.pdf",.plot_tft_1,height = 7, width = 5)

total_frequency_ns_by_cancertype_1kg = cancer_white %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type_ns(.,tally_1kg_white))) %>>%unnest()
.plot_bb_1 = bubble_plot_ns(total_frequency_ns_by_cancertype_1kg,total_frequency_ns_1kg)
ggsave("burden_plot/1kg/ns_bubble_plot_white.pdf",.plot_bb_1,height = 6,width = 8)

## UK10K
ns_focal_site_uk = quality_filter(norm_maf_all,.data_type="maf") %>>%
  filter(mutype == "missense"| mutype =="silent") %>>%
  dplyr::select(-alt) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),mutype=first(mutype)) %>>%
  ungroup() %>>%
  full_join(quality_filter(exac) %>>% filter(mutype == "missense"| mutype =="silent"))%>>%
  full_join(quality_filter(uk10k) %>>%filter(mutype == "missense"| mutype =="silent"))%>>%
  filter((is.na(AC_Adj)|AC_Adj/AN_Adj *100 <0.05), sum(AC_Adj,ac_cancer,ac_uk,na.rm = T) >2, chr!="chrX") %>>%
  mutate(focal = "ok") %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,focal)
uk_tft = uk10k %>>% quality_filter() %>>%
  dplyr::rename(ac_cont = ac_uk,an_cont = an_uk) %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cont,an_cont) %>>%
  filter(an_cont > max(an_cont)*0.5)
total_frequency_ns_uk =
  TFT_all_cancer_type_ns(cancer_white, uk_tft,.focal_site = ns_focal_site_uk)
.plot_bb_u = FDR_barplot_ns(total_frequency_ns_uk,.FDR = 0.2)
ggsave("burden_plot/uk10k/ns_fdr_white_likeCharles.pdf",.plot_tft_u,height = 7, width = 2)
total_frequency_ns_by_cancertype_uk = cancer_white %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,
                mutype,cancer_type,ac_cancer,an_cancer)%>>%
  nest(-cancer_type) %>>%
  mutate(data = purrr::map(data,~TFT_by_cancer_type_ns(.,uk_tft,ns_focal_site_uk))) %>>%
  unnest()
.plot_bb_u = bubble_plot_ns(total_frequency_ns_by_cancertype_uk,total_frequency_ns_uk)
ggsave("burden_plot/uk10k/ns_bubble_plot_white.pdf",.plot_bb_u,height = 6, width = 8)

cowplot::plot_grid(.plot_bb_w,.plot_bb_b,.plot_bb_1,.plot_bb_u,
                   ncol= 2)

####################################################################################################################
####################################################################################################################
####################################################################################################################
# nonsynonymou, synonymous
#white, black,1000genomes, UK10Kで比較
log_reverse = function(.log){
  parse(text = paste0("10^",-log10(.log)))
}
all_ns_fdr =  bind_rows(total_frequency_ns %>>%mutate(data_base = "ExAC_EUR"),
          total_frequency_ns_black %>>%mutate(data_base = "ExAC_AFR")) %>>%
  bind_rows(total_frequency_ns_1kg %>>%mutate(data_base = "1KG_EUR")) %>>%
  bind_rows(total_frequency_ns_uk  %>>%mutate(data_base = "UK10K")) 
all_ns_plot = function(.tbl, .mutype){
  .tbl$db_order = factor(.tbl$data_base,
                         levels = c("ExAC_EUR","ExAC_AFR","1KG_EUR","UK10K"))
  .plot = .tbl %>>%
    filter(mutype == .mutype) %>>%
    group_by(gene_symbol) %>>%
    filter(min(FDR) < 0.05) %>>% ungroup() %>>%
    mutate(FDR=10^(-log10(FDR)))%>>%
    mutate(FDR = ifelse(FDR>1e10,1e10,FDR)) %>>%
    ggplot(aes(x=gene_symbol, y=FDR ,fill=gene_symbol))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 20,size=0.2)+
    geom_hline(yintercept = 0,size=1,colour="gray")+
    scale_y_log10(labels = log_reverse, limits = c(1,1e10),expand = c(0,0))+
    #scale_fill_brewer(palette="Set1")+
    guides(fill="none")+
    facet_grid( db_order ~ mutype)+
    theme_bw()+
    theme( panel.grid.major.x = element_blank(),
           panel.background = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
  plot(.plot)
  return(.plot)
}
.n_all_fdr_plot = all_ns_fdr %>>%all_ns_plot(.mutype = "missense")
.s_all_fdr_plot = all_ns_fdr %>>%all_ns_plot(.mutype = "silent")
cowplot::plot_grid(.n_all_fdr_plot,.s_all_fdr_plot,ncol = 2,nrow = 1)
ggsave("burden_plot/fig/ns_fdr_all_db.pdf",height=12,width=8,
  cowplot::plot_grid(.n_all_fdr_plot,.s_all_fdr_plot+theme(axis.title = element_blank())
                     ,ncol = 2,nrow = 1))

#bubble plotも並べる
all_ns_fdr_by_cancertype =
  bind_rows(total_frequency_ns_by_cancertype,
            total_frequency_ns %>>%mutate(cancer_type ="Pan-cancer"))%>>%
  mutate(data_base = "ExAC_EUR") %>>%
  bind_rows(bind_rows(total_frequency_ns_by_cancertype_black,
                    total_frequency_ns_black %>>%mutate(cancer_type ="Pan-cancer"))%>>%
             mutate(data_base = "ExAC_AFR")) %>>%
  bind_rows(bind_rows(total_frequency_ns_by_cancertype_1kg,
                      total_frequency_ns_1kg %>>%mutate(cancer_type ="Pan-cancer")) %>>%
              mutate(data_base = "1KG_EUR")) %>>%
  bind_rows(bind_rows(total_frequency_ns_by_cancertype_uk,
             total_frequency_ns_uk %>>%mutate(cancer_type ="Pan-cancer")) %>>%
              mutate(data_base = "UK10K")) 
all_ns_bubble_plot = function(.tbl, .mutype){
  .tbl$db_order = factor(.tbl$data_base,
                         levels = c("ExAC_EUR","ExAC_AFR","1KG_EUR","UK10K"))
  .tbl = .tbl %>>%filter(mutype == .mutype)
  .plot = .tbl %>>%
    group_by(gene_symbol,data_base) %>>%
    filter(min(FDR) < 0.1) %>>% ungroup() %>>%
    filter(taf_cancer != 0)%>>%
    mutate(canty_order=ifelse(cancer_type=="Pan-cancer",2,1)) %>>%
    ggplot()+
    geom_point(aes(x=reorder(as.factor(cancer_type),canty_order),
                   y=reorder(as.factor(gene_symbol),desc(as.factor(gene_symbol))),
                   size=taf_cancer*100,color=cancer_type))+
    geom_point(data = .tbl%>>%filter(FDR < 0.05),
               aes(x=cancer_type,y=gene_symbol,size=taf_cancer*300),shape =21,stroke = 1)+
    #scale_color_brewer(palette = "Set1")+
    scale_size(breaks = c(10,20))+
    guides(colour="none",size= guide_legend("% of samples "))+
    facet_grid( db_order ~ mutype,scales = "free")+
    theme(panel.grid.minor = element_blank(),
          panel.background =  element_rect(fill = "white", colour = "grey50"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))
  plot(.plot)
  return(.plot)
}
.n_all_bubble_plot = all_ns_fdr_by_cancertype %>>%all_ns_bubble_plot(.mutype = "missense")
.s_all_bubble_plot = all_ns_fdr_by_cancertype %>>%all_ns_bubble_plot(.mutype = "silent")
cowplot::plot_grid(.n_all_bubble_plot,.s_all_bubble_plot,ncol = 2,nrow = 1)
ggsave("burden_plot/fig/ns_bubble_plot_all_db.pdf",height=12,width=12,
       cowplot::plot_grid(.n_all_bubble_plot+theme(legend.position = "none"),
                          .s_all_bubble_plot+theme(legend.position = "right")
                          ,ncol = 2,nrow = 1,rel_widths = c(1,1.4)))
  
  
  
  
  
#上同様にtotal frequency test
TFT_perm_by_row = function(.tbl){
  TFT = function(case_f,control_f,case,control){
    case_f=sum(case_f,na.rm = T)
    control_f=sum(control_f,na.rm = T)
    prop.test(c(case_f*case,control_f*control),c(case,control),alternative = "greater")$p.value
  }
  TFT_perm = function(maf_can,ac_cont,an_cont,perm=10000){
    #n_can=12896
    n_can=13018
    maf=sum(maf_can,na.rm = T)
    ac_cont=ac_cont[!is.na(ac_cont)]
    an_cont=an_cont[!is.na(an_cont)]
    maf_perm=rep(0,perm)
    for(i in 1:perm){
      for(g in 1:length(ac_cont)){
        maf_now=sum(sample.int(an_cont[g],n_can,replace=T)<ac_cont[g])/n_can
        maf_perm[i]=sum(maf_perm[i],maf_now)
      }
    }
    length(maf_perm[maf_perm>maf])/perm
  }
  #n_can=12896
  n_can=13018
  n_1kg=5008
  n_uk=4090
  n_exac=106210
  
  to_1kg = ifelse((sum(.tbl$maf_can)==0 & sum(.tbl$maf_1kg )==0),1,
                  ifelse((sum(.tbl$maf_can)>0&sum(.tbl$maf_1kg )>0),
                         TFT_perm(.tbl$maf_can,.tbl$ac_1kg  ,.tbl$an_1kg  ),
                         TFT(.tbl$maf_can,.tbl$maf_1kg ,n_can,n_1kg )))
  to_uk  = ifelse((sum(.tbl$maf_can)==0 & sum(.tbl$maf_uk  )==0),1,
                  ifelse((sum(.tbl$maf_can)>0&sum(.tbl$maf_uk  )>0),
                         TFT_perm(.tbl$maf_can,.tbl$ac_uk,.tbl$an_uk),
                         TFT(.tbl$maf_can,.tbl$maf_uk  ,n_can,n_uk  )))
  to_exac= ifelse((sum(.tbl$maf_can)==0 & sum(.tbl$maf_exac)==0),1,
                  ifelse((sum(.tbl$maf_can)>0&sum(.tbl$maf_exac)>0),
                         TFT_perm(.tbl$maf_can,.tbl$ac_exac ,.tbl$an_exac ),
                         TFT(.tbl$maf_can,.tbl$maf_exac,n_can,n_exac)))
  data.frame(to_1kg=to_1kg,to_uk=to_uk,to_exac=to_exac)
}

fortft_missense_silent = tally_norm_maf %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk,an_uk))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter(chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  mutate(maf_can  = ifelse(is.na(an_cancer),0,ifelse(ac_cancer/an_cancer>0.5,1-ac_cancer/an_cancer,ac_cancer/an_cancer)),
         maf_1kg  = ifelse(is.na(an_1kg   ),0,ifelse(ac_1kg   /an_1kg   >0.5,1-ac_1kg   /an_1kg   ,ac_1kg   /an_1kg   )),
         maf_uk   = ifelse(is.na(an_uk ),0,ifelse(ac_uk /an_uk >0.5,1-ac_uk /an_uk ,ac_uk /an_uk )),
         maf_exac = ifelse(is.na(an_exac  ),0,ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.05,mutype!="flank",mutype!="splice_region",mutype!="inframe_indel")%>>%
  mutate(MAF=class_separate(maf_exac))%>>%
  mutate(MAF_order=ifelse(MAF=="0~0.05",1,ifelse(MAF=="0.05~0.1",2,
                                                 ifelse(MAF=="0.1~0.5",3,ifelse(MAF=="0.5~1",4,5)))))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="missense"|mutype=="silent")
if(0){
  perm_tft0_5 = fortft_missense_silent %>>%
    filter(MAF_order<4)%>>% #MAF_order <4で0.5%以下,<3で0.1%以下,<2で0.05%以下
    nest(-gene_symbol,-mutype)%>>%
    mutate(tft=purrr::map(data,~TFT_perm_by_row(.))) %>>%
    dplyr::select(-data)%>>%unnest() %>>%
    mutate(`1000genomes`= p.adjust(to_1kg ,"fdr"),
           UK10K= p.adjust(to_uk  ,"fdr"),
           ExAC= p.adjust(to_exac,"fdr"))
  perm_tft0_05 = fortft_missense_silent %>>%
    filter(MAF_order<3)%>>% #MAF_order <4で0.5%以下,<3で0.1%以下,<2で0.05%以下
    nest(-gene_symbol,-mutype)%>>%
    mutate(tft=purrr::map(data,~TFT_perm_by_row(.))) %>>%
    dplyr::select(-data)%>>%unnest() %>>%
    mutate(`1000genomes`= p.adjust(to_1kg ,"fdr"),
           UK10K= p.adjust(to_uk  ,"fdr"),
           ExAC= p.adjust(to_exac,"fdr"))
}

chisq_tft0_5 = fortft_missense_silent %>>%
  filter(MAF_order<4)%>>% #MAF_order <4で0.5%以下,<3で0.1%以下,<2で0.05%以下
  group_by(gene_symbol,mutype,role) %>>%
  summarise(maf_can=sum(maf_can),maf_1kg=sum(maf_1kg),maf_uk=sum(maf_uk),maf_exac=sum(maf_exac)) %>>%
  nest(-gene_symbol,-mutype)%>>%
  mutate(tft=purrr::map(data,~TFT_by_row(.)))%>>%
  unnest()%>>%
  mutate(`1000genomes`= p.adjust(to_1kg ,"fdr"),
         UK10K= p.adjust(to_uk  ,"fdr"),
         ExAC= p.adjust(to_exac,"fdr"))
chisq_tft0_05 = fortft_missense_silent %>>%
  filter(MAF_order<2)%>>% #MAF_order <4で0.5%以下,<3で0.1%以下,<2で0.05%以下
  group_by(gene_symbol,mutype,role) %>>%
  summarise(maf_can=sum(maf_can),maf_1kg=sum(maf_1kg),maf_uk=sum(maf_uk),maf_exac=sum(maf_exac)) %>>%
  nest(-gene_symbol,-mutype)%>>%
  mutate(tft=purrr::map(data,~TFT_by_row(.)))%>>%
  unnest()%>>%
  mutate(`1000genomes`= p.adjust(to_1kg ,"fdr"),
         UK10K= p.adjust(to_uk  ,"fdr"),
         ExAC= p.adjust(to_exac,"fdr"))

.plot = chisq_tft0_05 %>>%
  filter(`1000genomes`<0.1 | UK10K<0.1|ExAC<0.1)%>>%
  tidyr::gather(database,FDR,`1000genomes`,UK10K,ExAC) %>>%
  mutate(FDR=FDR*(100^(-log10(FDR))))%>>%
  ggplot(aes(x=gene_symbol,y=FDR,fill=gene_symbol))+
  geom_bar(stat = "identity")+
  scale_y_log10(labels = log_reverse)+
  facet_grid(database ~ mutype,scales = "free")+
  guides(fill="none")+
  geom_hline(yintercept = 50,size=0.2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10,vjust = 0.5),axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12), axis.title.y = element_text(size=20),
        axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill="transparent", colour = "black"))
.plot
ggsave("burden_plot/silent_missense_tft_chisq_0.5.pdf",.plot,height = 7,width = 16)
ggsave("burden_plot/silent_missense_tft_chisq_0.05.pdf",.plot,height = 7,width = 12)

#作図(figure of allele frequency sum)
class_separate=function(.vec){
  ifelse(.vec < 0.0005, "0~0.05",
         ifelse(.vec < 0.001,"0.05~0.1",
                ifelse(.vec < 0.005,"0.1~0.5",
                       ifelse(.vec <0.01,"0.5~1","1~5"))))
}

site_number = tally_norm_maf %>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk,an_uk))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter(chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  mutate(maf_exac = ifelse(is.na(an_exac  ),0,ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.05,mutype!="flank",mutype!="splice_region",mutype!="inframe_indel")%>>%
  mutate(MAF=class_separate(maf_exac),mutype=ifelse(mutype=="splice", "truncating", mutype))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene"))%>>%
  group_by(role,mutype,MAF)%>>%
  summarise(cancer=sum(!is.na(an_cancer)),`1000genomes`=sum(!is.na(an_1kg)),
            UK10K=sum(!is.na(an_uk)),ExAC=sum(!is.na(an_exac))) %>>%
  tidyr::gather(data_base,site_num,cancer,`1000genomes`,UK10K,ExAC)

.plot = tally_norm_maf %>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk,an_uk))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter(chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  mutate(maf_can  = ifelse(is.na(an_cancer),0,
                           ifelse(ac_cancer/an_cancer>0.5,1-ac_cancer/an_cancer,ac_cancer/an_cancer)),
         maf_1kg  = ifelse(is.na(an_1kg   ),0,
                           ifelse(ac_1kg   /an_1kg   >0.5,1-ac_1kg   /an_1kg   ,ac_1kg   /an_1kg   )),
         maf_uk   = ifelse(is.na(an_uk ),0,
                           ifelse(ac_uk /an_uk >0.5,1-ac_uk /an_uk ,ac_uk /an_uk )),
         maf_exac = ifelse(is.na(an_exac  ),0,
                           ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.05,mutype!="flank",mutype!="splice_region",mutype!="inframe_indel")%>>%
  mutate(MAF=class_separate(maf_exac),mutype=ifelse(mutype=="splice", "truncating", mutype))%>>%
  mutate(MAF_order=ifelse(MAF=="0~0.05",1,ifelse(MAF=="0.05~0.1",2,
                                                 ifelse(MAF=="0.1~0.5",3,ifelse(MAF=="0.5~1",4,5)))))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene"))%>>%
  group_by(role,mutype,MAF,MAF_order)%>>%
  summarise(cancer=sum(maf_can),`1000genomes`=sum(maf_1kg),
            UK10K=sum(maf_uk),ExAC=sum(maf_exac)) %>>%
  tidyr::gather(data_base,sum_of_AF,cancer,`1000genomes`,UK10K,ExAC)%>>%
  mutate(data_base_order=ifelse(data_base=="cancer",1,
                                ifelse(data_base=="1000genomes",2,
                                       ifelse(data_base=="UK10K",4,3))))%>>%
  left_join(site_number)%>>%
  ggplot(aes(x=reorder(MAF,MAF_order),y=sum_of_AF,fill=reorder(data_base,data_base_order)))+
  geom_bar(stat ="identity",position = "dodge")+
  facet_grid(mutype ~ role,scales = "free")+
  labs(x="MAF (%)", y="sum of MAF (%)", fill="data base")+
  geom_text(aes(group=reorder(data_base,data_base_order),hjust=-0.1*log10(site_num),label=site_num),
            position = position_dodge(width = 0.9), size = 3,angle=90)+
  theme_bw()+
  theme(axis.title = element_text(size=30), axis.text.x = element_text(angle=-45, size=20,hjust = 0),
        axis.text.y = element_text(size = 20),strip.text = element_text(size = 20),
        panel.grid.major.x = element_blank(), strip.background = element_rect(fill="transparent", colour = "black"))
.plot
ggsave("burden_plot/sum_AF_barplot.pdf",.plot,height = 15,width = 12)


