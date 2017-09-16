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

tally_1kg_all= maf_1kg_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_1kg=n()) %>>% ungroup() %>>%
  left_join(maf_1kg_all %>>% group_by(sample,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_1kg=n())) %>>%
  mutate(an_1kg = ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg))%>>%
  mutate(chr=paste0("chr",chr),homo_1kg=ifelse(is.na(homo_1kg),0,homo_1kg)) #%>>%
  #full_join(vcf_1kg_ac0)

rm(sample_num_1kg,sample_num_X_1kg)
####UK 10K ###
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle_likevcf.tsv.gz",col_types = "cdccdddd")
maf_10k=strip_maf("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle.maf")
uk10k=maf_10k %>>%
  left_join(vcf_10k%>>%mutate(start=ifelse(ref=="-",start -1,start))) %>>%
  mutate(chr=paste0("chr",chr)) %>>%
  dplyr::rename(ac_uk10k=uk_ac,an_uk10k=uk_an)%>>%
  filter(an_uk10k >3000)
rm(vcf_10k,maf_10k)

### EXAC nonTCGA ###
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",
                  col_types = "cdccdd") %>>%dplyr::rename(start=posi) 
exac_nonindel=strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver.maf") %>>%
  left_join(vcf_exac)
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                        col_types = "cdccdd") %>>%dplyr::rename(start=posi)
exac_indel = strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver_indel.maf") %>>%
  left_join(vcf_exac_indel)
vcf_exac =rbind(vcf_exac,vcf_exac_indel) %>>%mutate(chr=paste0("chr",chr))
exac = rbind(exac_nonindel,exac_indel) %>>%mutate(chr=paste0("chr",chr))
rm(vcf_exac_indel,exac_nonindel,exac_indel)


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#ここから解析
#kokokarakaiseki
#here start analysis
#####################################################################################################################
### Charlesら(NatCom2015の解析の確認)
TFT_by_row = function(.tbl){
  TFT = function(case_f,control_f,case,control){
    case_m = round(case * case_f)
    control_m = round(control * control_f)
    prop.test(c(case_m,control_m),c(case,control),alternative = "greater")$p.value
  }
  n_can=12704
  n_1kg=5008
  n_uk=4090
  n_exac=106210
  .tbl %>>%
    mutate(to_1kg = ifelse((maf_can==0 & maf_1kg ==0),1, TFT(maf_can,maf_1kg ,n_can,n_1kg )),
           to_uk  = ifelse((maf_can==0 & maf_uk  ==0),1, TFT(maf_can,maf_uk  ,n_can,n_uk  )),
           to_exac= ifelse((maf_can==0 & maf_exac==0),1, TFT(maf_can,maf_exac,n_can,n_exac)))%>>%
    dplyr::select(to_1kg,to_uk,to_exac)
}

total_frequency_truncate = tally_norm_maf %>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk10k,an_uk10k))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter((mutype=="truncating"|mutype=="splice") & chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  #mutate(ac_cancer= ifelse(is.na(ac_cancer),0,ac_cancer),
  #       ac_1kg  = ifelse(is.na(ac_1kg   ),0,ac_1kg   ),
  #       ac_uk10k= ifelse(is.na(ac_uk10k ),0,ac_uk10k ),
  #       ac_exac = ifelse(is.na(ac_exac  ),0,ac_exac  )) %>>%
  #mutate(ac_all=ac_cancer + ac_1kg + ac_uk10k + ac_exac) %>>%
  #filter(ac_all >3) %>>% #ここをいれるとCharlesらと同じになる
  #filter(!(ac_cancer<2 & ac_1kg<2 & ac_uk10k<2 & ac_exac <2)) %>>%
  mutate(maf_can  = ifelse(is.na(an_cancer),0,ifelse(ac_cancer/an_cancer>0.5,1-ac_cancer/an_cancer,ac_cancer/an_cancer)),
         maf_1kg  = ifelse(is.na(an_1kg   ),0,ifelse(ac_1kg   /an_1kg   >0.5,1-ac_1kg   /an_1kg   ,ac_1kg   /an_1kg   )),
         maf_uk   = ifelse(is.na(an_uk10k ),0,ifelse(ac_uk10k /an_uk10k >0.5,1-ac_uk10k /an_uk10k ,ac_uk10k /an_uk10k )),
         maf_exac = ifelse(is.na(an_exac  ),0,ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.0005)%>>%
  group_by(gene_symbol)%>>%
  summarise(maf_can=sum(maf_can),maf_1kg=sum(maf_1kg),maf_uk=sum(maf_uk),maf_exac=sum(maf_exac))%>>%
  nest(-gene_symbol) %>>%
  mutate(tft=purrr::map(data,~TFT_by_row(.)))%>>%
  unnest()%>>%
  mutate(`1000genomes`= p.adjust(to_1kg ,"fdr"),
         UK10K= p.adjust(to_uk  ,"fdr"),
         ExAC= p.adjust(to_exac,"fdr"))

log_reverse=function(l){
  l <-format(l,scientific = T)
  l <- gsub("\\+","-",l)
  l <- gsub("e", "0^", l)
}
.plot = total_frequency_truncate %>>%
  filter(`1000genomes`<0.1 | UK10K<0.1|ExAC<0.1)%>>%
  tidyr::gather(database,FDR,`1000genomes`,UK10K,ExAC) %>>%
  mutate(FDR=FDR*(100^(-log10(FDR))))%>>%
  ggplot(aes(x=gene_symbol,y=FDR))+
  geom_bar(stat = "identity")+
  scale_y_log10(labels = log_reverse)+
  facet_grid(database ~ .,scales = "free")
ggsave("burden_plot/truncate_fdr_likeCharles.pdf",.plot,height = 10,width = 10)

####### ↑を確認用にCAST法でやってみる #############

####################################################################################################################
#作図(figure of allele frequency sum)
class_separate=function(.vec){
  ifelse(.vec < 0.0005, "0~0.05%",
         ifelse(.vec < 0.001,"0.05~0.1%",
                ifelse(.vec < 0.005,"0.1~0.5%",
                       ifelse(.vec <0.01,"0.5~1%","1~5%"))))
}

.plot = tally_norm_maf %>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk10k,an_uk10k))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter(chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  mutate(maf_can  = ifelse(is.na(an_cancer),0,ifelse(ac_cancer/an_cancer>0.5,1-ac_cancer/an_cancer,ac_cancer/an_cancer)),
         maf_1kg  = ifelse(is.na(an_1kg   ),0,ifelse(ac_1kg   /an_1kg   >0.5,1-ac_1kg   /an_1kg   ,ac_1kg   /an_1kg   )),
         maf_uk   = ifelse(is.na(an_uk10k ),0,ifelse(ac_uk10k /an_uk10k >0.5,1-ac_uk10k /an_uk10k ,ac_uk10k /an_uk10k )),
         maf_exac = ifelse(is.na(an_exac  ),0,ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.05,mutype!="flank",mutype!="splice_region",mutype!="inframe_indel")%>>%
  mutate(MAF=class_separate(maf_exac))%>>%
  mutate(MAF_order=ifelse(MAF=="0~0.05%",1,ifelse(MAF=="0.05~0.1%",2,ifelse(MAF=="0.1~0.5%",3,ifelse(MAF=="0.5~1%",4,5)))))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene"))%>>%
  group_by(role,mutype,MAF,MAF_order)%>>%
  summarise(cancer=sum(maf_can),`1000genomes`=sum(maf_1kg),
            UK10K=sum(maf_uk),ExAC=sum(maf_exac)) %>>%
  tidyr::gather(data_base,sum_of_AF,cancer,`1000genomes`,UK10K,ExAC)%>>%
  mutate(data_base_order=ifelse(data_base=="cancer",1,ifelse(data_base=="1000genomes",2,
                                                             ifelse(data_base=="UK10K",3,4))))%>>%
  ggplot(aes(x=reorder(MAF,MAF_order),y=sum_of_AF,fill=reorder(data_base,data_base_order)))+
  geom_bar(stat ="identity",position = "dodge")+
  facet_grid(mutype ~ role,scales = "free")
ggsave("burden_plot/sum_AF_barplot.pdf",.plot,height = 15,width = 15)

#上同様にtotal frequency test
fortft_missense_silent = tally_norm_maf %>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk10k,an_uk10k))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter(chr!="chrX")%>>%
  remove_duplicate()%>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n))%>>%dplyr::select(-n)%>>%
  mutate(maf_can  = ifelse(is.na(an_cancer),0,ifelse(ac_cancer/an_cancer>0.5,1-ac_cancer/an_cancer,ac_cancer/an_cancer)),
         maf_1kg  = ifelse(is.na(an_1kg   ),0,ifelse(ac_1kg   /an_1kg   >0.5,1-ac_1kg   /an_1kg   ,ac_1kg   /an_1kg   )),
         maf_uk   = ifelse(is.na(an_uk10k ),0,ifelse(ac_uk10k /an_uk10k >0.5,1-ac_uk10k /an_uk10k ,ac_uk10k /an_uk10k )),
         maf_exac = ifelse(is.na(an_exac  ),0,ifelse(ac_exac  /an_exac  >0.5,1-ac_exac  /an_exac  ,ac_exac  /an_exac  ))) %>>%
  filter(maf_exac<0.05,mutype!="flank",mutype!="splice_region",mutype!="inframe_indel")%>>%
  mutate(MAF=class_separate(maf_exac))%>>%
  mutate(MAF_order=ifelse(MAF=="0~0.05%",1,ifelse(MAF=="0.05~0.1%",2,
                    ifelse(MAF=="0.1~0.5%",3,ifelse(MAF=="0.5~1%",4,5)))))%>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  group_by(gene_symbol,mutype,MAF,MAF_order)%>>%
  summarise(cancer=sum(maf_can),`1000genomes`=sum(maf_1kg),
            UK10K=sum(maf_uk),ExAC=sum(maf_exac)) %>>%
  filter(mutype=="missense"|mutype=="silent")
