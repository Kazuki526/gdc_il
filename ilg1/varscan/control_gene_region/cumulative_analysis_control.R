######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.R, quality_filter.Rを実行してからこのスクリプトを開始する！！！


######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
# ref_minor = ref_minor_focal_cont %>>%
#   filter(focal=="ok") %>>%dplyr::select(-focal)%>>%
#   left_join(left_join(tally_norm_maf_cont,vcf_exac_cont%>>%
#                         mutate(AF=AC_Adj/AN_Adj)%>>%
#                         dplyr::select(chr,start,ref,alt,AF))%>>%
#               dplyr::select(-ac_cancer,-hom_cancer)) %>>%
#   filter(AF > 0.5) %>>% mutate(MAF=1-AF) %>>%
#   left_join(norm_maf_all_cont) %>>%
#   mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
#          n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
#          soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
#          LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
#   dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
#                 mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,age,Protein_position) %>>%
#   quality_filter_cont(.data_type = "maf",.varscan = T) %>>%
#   mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
#   filter(MAC > 0) %>>%
#   dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age)
# 
# all_maf_for_cumulative_cont = norm_maf_all_cont %>>%
#   quality_filter_cont(.data_type = "maf",.varscan = T) %>>%
#   left_join(vcf_exac_cont%>>%mutate(AF=AC_Adj/AN_Adj)%>>%dplyr::select(chr,start,ref,alt,AF)) %>>%
#   mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
#   filter(AF < 0.5) %>>% mutate(MAF = AF) %>>%
#   mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
#   dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age) %>>%
#   rbind(ref_minor) %>>%
#   filter(!is.na(age)) %>>%
#   filter(!is.na(age),chr != "chrX")
#write_df(all_maf_for_cumulative_cont,"all_patient/all_maf_for_cumulative_control.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
#####################################################################################################################
#####################################################################################################################
#violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
truncating_filling = function(.truncating_count){
  .max_count=max(.truncating_count$truncating_count_n)
  .max_count=ifelse(.max_count>9,9,.max_count)
  .outcome = tibble::tibble(truncating_count_n=1:.max_count) %>>%
    left_join(.truncating_count%>>%count(truncating_count_n)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$truncating_count_n) != 0){
    .outcome = .outcome %>>%
      mutate(truncating_num = truncating_count_n,
             truncating_num_order = truncating_count_n,
        age=-50, cancer_type="BRCA", patient_id="no_patient",gender="male")
    return(rbind(.truncating_count,.outcome))
  }else{return(.truncating_count)}
}
#####################################################################################################################
#指数表記
exponent_notation = function(.num){
  .log=trunc(log10(.num))
  .log=ifelse(.log>0,.log+1,.log-1)
  paste0(.num*(10^-.log)," %*% 10^",.log)
}

#####################################################################################################################
#ここから解析
#ここから解析
#ここから解析

###################################################################################################################
####################################################### TSG #######################################################
###################################################################################################################
perm_pvalue = function(.regression,.CT=NA,.tbl,.tail){
  if(!is.na(.CT)){.tbl=filter(.tbl,cancer_type==.CT)}
  if(.tail == "one"){
    .tbl = .tbl %>>%
      filter(regression < .regression )
  }else if(.tail == "two"){
    .tbl = .tbl %>>%
      filter(regression > abs(.regression))
  }
  length(.tbl$regression) /10000
}
truncate_plot_allcantype= function(.tbl,.permu=T,.test_tail="one",.permu_file=NA){
  .max_count=ifelse(max(.tbl$truncating_count_n) >10,10,max(.tbl$truncating_count_n))
  .p_posi=.max_count/2 + 1
  .coef_posi=.max_count+1
  lm=lm(age ~ truncating_count_n, data=.tbl)
  regression=as.data.frame(as.list(coef(lm)))
  trunc_permu = function(.times,.tbl){
    if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
    .tbl_sample=.tbl%>>%mutate(age=sample(age,length(age)))
    as.tibble(as.list(coef(lm(age ~ truncating_count_n,data = .tbl_sample))))
  }
  if(.permu){
    if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
      regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
    }else{
      regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
        mutate(tbl=purrr::map(times,~trunc_permu(.times = .,.tbl = .tbl))) %>>%
        unnest() %>>%
        rename(regression = truncating_count_n)
      write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
    }
    .p=perm_pvalue(.tbl = regression_tbl,
                   .regression = regression$truncating_count_n,
                   .tail = .test_tail)
    regression = mutate(regression,p_value = .p)%>>%(?.)%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",p_value))
                             ,ifelse(p_value<0.05,"*","")))
  }else{
    regression = regression %>>%
      mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                             summary(lm)$fstatistic["dendf"]))) %>>%(?.)%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",p_value))
                             ,ifelse(p_value<0.05,"*","")))
  }
  .tbl = .tbl %>>%
    mutate(truncating_num=as.character(ifelse(truncating_count_n >= 10,"10-",truncating_count_n)),
           truncating_num_order=ifelse(truncating_count_n >=10,10,truncating_count_n)) %>>%
    truncating_filling() %>>%
    mutate(cancer_type="All Cancer Type")
  .plot = .tbl %>>%
    ggplot(aes(x=reorder(truncating_num,truncating_num_order), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    ylim(0,90)+
    geom_text(data =.tbl %>>%count(truncating_num),
              aes(x=truncating_num,y=5,label=n),size=4,position="stack")+
    geom_abline(data = regression %>>%filter(p_value > 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n),hjust=1)+
    geom_text(data = regression,aes(x=.coef_posi,y=10,label=out_reg),size=7,hjust=1)+
    facet_wrap( ~ cancer_type)+
    xlab("Number of Truncated Gene")+ylab("Age at Diagnosis")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}
###########################################################################################
truncate_plot_bycantype = function(.tbl,.permu=T,.test_tail="one",.permu_file=NA){
  .max_count=ifelse(max(.tbl$truncating_count_n) >10,10,max(.tbl$truncating_count_n))
  .p_posi=.max_count/2 + 1
  .coef_posi=.max_count+1
  lm_p=function(.data,.permu){
    if(max(.data$truncating_count_n) == 0){
      .dtbl = tibble::tibble(X.Intercept. =0, truncating_count_n =0)
      if(.permu){return(.dtbl)}else{return(mutate(.dtbl,p_value=1))}
    }else{
      lm=lm(age ~ truncating_count_n, data=.data)
      .dtbl = as.data.frame(as.list(coef(lm)))
      if(.permu){return(.dtbl)
      }else{
        return(mutate(.dtbl,p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                                            summary(lm)$fstatistic["dendf"]))))
      }
    }
  }
  trunc_permu_byCT = function(.times,.tbl){
    if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
    .tbl%>>%mutate(age=sample(age,length(age)))%>>%
      nest(-cancer_type)%>>%dplyr::rename(data_=data)%>>%
      mutate(data_ = purrr::map(data_,~lm_p(.,.permu = T)))%>>%unnest()
  }
  regression = .tbl %>>%
    tidyr::nest(-cancer_type) %>>%
    mutate(data = purrr::map(data,~lm_p(.,.permu = .permu))) %>>%
    unnest()
  if(.permu){#if do permutation
    if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
      regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
    }else{
      regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
        mutate(data = purrr::map(times,~trunc_permu_byCT(.times =.,.tbl = .tbl))) %>>%
        unnest() %>>%
        rename(regression = truncating_count_n)
      write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
    }
    regression = regression %>>%
      mutate(p_value = pmap_dbl(.,function(cancer_type,truncating_count_n,...){
        perm_pvalue(.regression =truncating_count_n,.CT=cancer_type,
                    .tbl = regression_tbl,.tail = .test_tail)})) %>>%(?.)
  }
  .tbl =.tbl %>>%
    mutate(truncating_num=as.character(ifelse(truncating_count_n >= 10,"10-",truncating_count_n)),
           truncating_num_order=ifelse(truncating_count_n >=10,10,truncating_count_n)) %>>%
    truncating_filling()
  .plot = .tbl %>>%
    ggplot(aes(x=reorder(truncating_num,truncating_num_order), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    ylim(0,90)+
    facet_wrap( ~ cancer_type, ncol = 5)+
    geom_text(data =.tbl %>>%count(cancer_type,truncating_num),
              aes(x=truncating_num,y=5,label=n),size=2.5,position="stack",angle=-45)+
    geom_abline(data = regression %>>%filter(p_value > 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),size=5,hjust=1)+
    geom_text(data = regression %>>%mutate(coefi=paste0("regression=",round(truncating_count_n*100)/100)),
              aes(x=.coef_posi,y=15,label =coefi),size=3.5,hjust=1)+
    xlab("number of truncated gene")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=10),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}
#患者ごとのtruncating な遺伝子の数
truncating_count_cont = all_maf_for_cumulative_cont %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)


.plot_all = truncating_count_cont %>>%
  truncate_plot_allcantype(.permu_file = "control/truncate_all.tsv")
ggsave("age_plot/control_region/truncating.pdf",.plot_all,height = 5,width = 8)
.plot_by = truncating_count_cont%>>%
  truncate_plot_bycantype(.permu_file = "control/truncate_all_byCT.tsv")
ggsave("age_plot/control_region/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p_value=0.01401 しかしtruncating=0の患者が少ないから、、
t.test(truncating_count_cont[truncating_count_cont$truncating_count_n>0,]$age/365.25,
       truncating_count_cont[truncating_count_cont$truncating_count_n==0,]$age/365.25,alternative="less")
####### figure 用
.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/control.pdf",.plot,width = 14,height = 8)


############### MAF<0.05%のtruncating mutationのみでやってみたら？
truncating_count_cont_rare = all_maf_for_cumulative_cont %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(MAF < 0.05)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)


.plot_all = truncating_count_cont_rare%>>%
  truncate_plot_allcantype(.permu_file = "control/truncate_rare.tsv")
ggsave("age_plot/control_region/rare_truncating.pdf",.plot_all,height = 5,width = 8)
.plot_by = truncating_count_cont_rare%>>%
  truncate_plot_bycantype(.permu_file = "control/truncate_rare_byCT.tsv")
ggsave("age_plot/control_region/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.2341
t.test(truncating_count_cont_rare[truncating_count_cont_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_cont_rare[truncating_count_cont_rare$truncating_count_n==0,]$age/365.25,alternative="less")
####### figure 用
.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/control_rare.pdf",.plot,width = 14,height = 8)
