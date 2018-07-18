all_maf_for_cumulative = read_tsv("all_patient/all_maf_for_cumulative.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
all_patient_info = read_tsv("~/git/all_patient/all_patient_response.tsv")
####################################################################################################################
####################################################################################################################
####################################################################################################################
all_maf_for_cumulative_cont %>>%
  group_by(gene_symbol,chr,ref,alt,start,end,MAF,mutype) %>>%
  summarise(MAC=sum(MAC)) %>>%
  filter(mutype=="silent",MAF>0.06) %>>%View
regression_plot_by_binsize = function(.maf,.role = NULL,.mutype="missense",.max_MAF = 1){
  regression_out = function(.maf){
    ##missense の数
    missense_count = .maf %>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>%
      {left_join(patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num))
    #相関直線を
    lm=lm(age/365.25 ~ missense_num, data=missense_count)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"]))
  }
  if(!is.null(.role)){
    .maf = .maf %>>%left_join(driver_genes%>>%dplyr::select(gene,role),by = c("gene_symbol"="gene"))%>>%
      filter(role == .role)
  }
  .maf = mutate(.maf,MAF=MAF*100)
  #１つのbinに何site入れるか(MAF==0のsite数にしている)
  .bin_mac_count = .maf %>>%
    filter(mutype==.mutype, MAF==0) %>>%
    {sum(.$MAC)}
  print(paste0("bin size mutatin num = ",.bin_mac_count))
  #上のbin sizeでグループ分け
  .site_group_by_MAF = .maf %>>%
    filter(MAF<.max_MAF, mutype == .mutype) %>>%
    group_by(MAF,gene_symbol,chr,start,end,ref,alt,mutype) %>>%
    summarise(MAC=sum(MAC)) %>>% ungroup() %>>%
    dplyr::arrange(MAF) %>>%
    mutate(cumsum_MAC = cumsum(MAC)) %>>%
    mutate(bin_group = cumsum_MAC %/% .bin_mac_count) %>>%
    dplyr::select(-MAC,cumsum_MAC) %>>%
    group_by(bin_group) %>>%
    mutate(max_MAF = max(MAF), min_MAF = min(MAF)) %>>%
    ungroup()
  .first_bin = first(.site_group_by_MAF$max_MAF[.site_group_by_MAF$bin_group==1])/2#plot作製用
  #plot
  .for_plot = .maf %>>%
    filter(mutype == .mutype) %>>%
    left_join(.site_group_by_MAF) %>>%
    tidyr::nest(-bin_group,-max_MAF,-min_MAF) %>>%
    mutate(tbl = purrr::map(data, ~regression_out(.)))%>>%
    dplyr::select(-data) %>>%
    unnest() %>>%
    mutate(min_MAF = ifelse(bin_group ==0,-.first_bin,min_MAF),
           significance = ifelse(p_value < 0.05,"significance","not significance"))
  .plot = ggplot()+
    geom_rect(aes(xmin=min_MAF,xmax=max_MAF, ymin=0,ymax=missense_num,fill = significance),
              data = .for_plot %>>%filter(missense_num>0))+
    geom_rect(aes(xmin=min_MAF,xmax=max_MAF, ymax=0,ymin=missense_num,fill = significance),
              data = .for_plot %>>%filter(missense_num<0))+
    scale_fill_manual(values = c("grey","red") )
  plot(.plot)
  return(.plot)
}

if(exists("all_maf_for_cumulative")){
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "TSG",.mutype = "missense") 
  ggsave("age_plot/regression_by_MAF_group/TSG_missense_1.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "TSG",
                                     .mutype = "missense", .max_MAF = 10) 
  ggsave("age_plot/regression_by_MAF_group/TSG_missense_10.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "TSG",.mutype = "silent")
  ggsave("age_plot/regression_by_MAF_group/TSG_silent_1.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "TSG",
                                     .mutype = "silent",.max_MAF = 10)
  ggsave("age_plot/regression_by_MAF_group/TSG_silent_10.pdf",.plot)
  
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "oncogene",.mutype = "missense") 
  ggsave("age_plot/regression_by_MAF_group/oncogene_missense_1.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "oncogene",
                                     .mutype = "missense", .max_MAF = 10) 
  ggsave("age_plot/regression_by_MAF_group/oncogene_missense_10.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "oncogene",.mutype = "silent")
  ggsave("age_plot/regression_by_MAF_group/oncogene_silent_1.pdf",.plot)
  .plot = regression_plot_by_binsize(all_maf_for_cumulative,.role = "oncogene",
                                     .mutype = "silent",.max_MAF = 10)
  ggsave("age_plot/regression_by_MAF_group/oncogene_silent_10.pdf",.plot)
}


if(exists("all_maf_for_cumulative_cont")){
  .plot = regression_plot_by_binsize(all_maf_for_cumulative_cont)
  ggsave("age_plot/regression_by_MAF_group/contol_missense_1.pdf")
  .plot = regression_plot_by_binsize(all_maf_for_cumulative_cont,.max_MAF = 10)
  ggsave("age_plot/regression_by_MAF_group/contol_missense_10.pdf")
  .plot = regression_plot_by_binsize(all_maf_for_cumulative_cont,.mutype = "silent")
  ggsave("age_plot/regression_by_MAF_group/contol_silent_1.pdf")
  .plot = regression_plot_by_binsize(all_maf_for_cumulative_cont,.mutype = "silent",.max_MAF = 10)
  ggsave("age_plot/regression_by_MAF_group/contol_silent_10.pdf")
}
