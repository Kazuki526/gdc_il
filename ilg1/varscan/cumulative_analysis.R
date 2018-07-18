######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.R, quality_filter.Rを実行してからこのスクリプトを開始する！！！
loadNamespace('cowplot')

######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
maf_trim_for_cumulative = function(.maf=norm_maf_all,.vcf=vcf_exac,.race="all",.fdr=0.01,
                                   .database="all",.duplicate=T,.somatic=T,.varscan=T){
  if(.race=="white"){
    .maf = .maf %>>%left_join(patient_race) %>>%filter(race==.race)
    .vcf = .vcf %>>%mutate(AF=(AC_FIN+AC_NFE)/(AN_FIN+AN_NFE)) %>>%
      dplyr::select(chr,start,ref,alt,AF)
  }else if(.race=="black"){
    .maf = .maf %>>%left_join(patient_race) %>>%filter(race==.race)
    .vcf = .vcf %>>%mutate(AF=(AC_AFR)/(AN_AFR)) %>>%
      dplyr::select(chr,start,ref,alt,AF) 
  }else if(.race=="all"){
    if(substitute(.vcf)=="vcf_exac"){
      .vcf = .vcf %>>%
        mutate(AF=AC_Adj/AN_Adj)%>>%
        dplyr::select(chr,start,ref,alt,AF)
    }}else {stop(paste0(".race is wrong .race=",.race,"\nuse all, white or black!"))}
  ref_minor = mid_af_coverage %>>%
    left_join(patient_race) %>>%
    filter(if(.race !="all"){race==.race}else{chr==chr}) %>>%
    left_join(left_join(tally_norm_maf,.vcf)%>>%
                dplyr::select(-ac_cancer,-hom_cancer)) %>>%
    filter(AF > 0.5) %>>% mutate(MAF=1-AF) %>>%
    left_join(.maf) %>>%
    mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
           n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
           soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
           LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
                  mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,age,Protein_position) %>>%
    quality_filter(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
    filter(MAC > 0) %>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age,Protein_position)
  
  .maf %>>%
    quality_filter(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    left_join(.vcf) %>>%
    mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
    filter(AF < 0.5) %>>% mutate(MAF = AF) %>>%
    mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age,Protein_position) %>>%
    rbind(ref_minor) %>>%
    filter(!is.na(age))
}

#all_maf_for_cumulative = maf_trim_for_cumulative()
#write_df(all_maf_for_cumulative,"all_patient/all_maf_for_cumulative.tsv.gz")
all_maf_for_cumulative = read_tsv("all_patient/all_maf_for_cumulative.tsv.gz")

#####################################################################################################################
#####################################################################################################################
#violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
truncating_filling = function(.truncating_count){
  .max_count=max(.truncating_count$truncating_count_n)
  .outcome = data.frame(truncating_count_n=1:.max_count) %>>%
    left_join(.truncating_count%>>%count(truncating_count_n)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$truncating_count_n) != 0){
    .outcome = .outcome %>>%
      mutate(age=-50, cancer_type="BRCA", patient_id="no_patient",gender="male")
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

####################################################################################################################
####################################################### TSG ########################################################
####################################################################################################################
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

truncate_plot_allcantype= function(.tbl,.permu=T,.test_tail="one",.permu_do=F,.permu_file=NA){
  .max_count=max(.tbl$truncating_count_n)
  .coef_posi=ifelse(.max_count==1,2.25,ifelse(.max_count==2,3,.max_count))
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
                             ifelse(p_value==0,"<0.0001",paste0("=",p_value))))
  }else{
    regression = regression %>>%
      mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                             summary(lm)$fstatistic["dendf"]))) %>>%(?.)%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  }
  .plot = .tbl %>>%truncating_filling()%>>%
    mutate(cancer_type ="All Cancer Types") %>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    geom_text(data =.tbl %>>%count(truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=5,label=n),size=3,position="stack")+
    geom_abline(data = regression %>>%filter(p_value > 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n))+
    geom_text(data = regression, aes(x=.coef_posi,y=10,label=out_reg),size=5,hjust=1)+
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
##############################################################################################
truncate_plot_bycantype = function(.tbl,.permu=T,.test_tail="one",.permu_do=F,.permu_file=NA){
  .max_count=max(.tbl$truncating_count_n)
  .coef_posi=ifelse(.max_count == 1,2.4,.max_count+1)
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
                    .tbl = regression_tbl,.tail = .test_tail)})) %>>%(?.)%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  }
  .plot = .tbl %>>%truncating_filling()%>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    facet_wrap( ~ cancer_type, ncol = 5)+
    geom_text(data =.tbl %>>%count(cancer_type,truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=2.5,label=n),size=2.5,position="stack",angle=-45)+
    geom_abline(data = regression %>>%filter(p_value > 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n))+
    geom_text(data = regression, aes(x=.coef_posi,y=15,label =out_reg),size=3.5,hjust=1)+
    xlab("Number of Truncated Gene")+ylab("Age at Diagnosis")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}

#####################################################################################################
#患者ごとのtruncating な遺伝子の数
truncating_count = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%#filter(truncating_count_n<5)%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all.tsv")
ggsave("age_plot/cumulative/truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT.tsv")
ggsave("age_plot/cumulative/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p_value=0.1965
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
truncating_count_rare = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

options(scipen=1)
.plot_all = truncating_count_rare%>>%#filter(truncating_count_n<5)%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_rare.tsv") 
ggsave("age_plot/cumulative/rare_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_rare%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.02574
t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/TSG_rare.pdf",.plot,width = 14,height = 8)

if(0){
##### rare BRCA1,2を別にしてみたら？####
brca_truncate_num = patient_list %>>%
  left_join(all_maf_for_cumulative %>>%
              filter((gene_symbol=="BRCA1"|gene_symbol=="BRCA2")&MAF<0.0005) %>>%  ####ここを変えてBRCA1,2のみか以外か変える
              filter((mutype=="truncating"|mutype=="splice")) %>>%
              count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
              group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot_all = brca_truncate_num %>>%truncate_plot_allcantype()
ggsave("age_plot/cumulative/brca_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = brca_truncate_num %>>%truncate_plot_bycantype()
ggsave("age_plot/cumulative/brca_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/rare_brca12.pdf",.plot,width = 10,height = 10)

##filterをいじってBRCA1,2以外では？
brca_truncate_num_ = patient_list %>>%
  left_join(all_maf_for_cumulative %>>%
              filter(!((gene_symbol=="BRCA1"|gene_symbol=="BRCA2")&MAF<0.0005)) %>>%  ####ここを変えてBRCA1,2のみか以外か変える
              left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
              filter(mutype=="truncating"|mutype=="splice") %>>%
              filter(role=="TSG"|role=="oncogene/TSG")%>>%
              count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
              group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)
.plot = brca_truncate_num_ %>>%truncate_plot_allcantype()
ggsave("age_plot/cumulative/brca_non_truncating.pdf",.plot,height = 5,width = 5)
.plot = brca_truncate_num_ %>>%truncate_plot_bycantype()
ggsave("age_plot/cumulative/brca_non_truncating_by_cancerype.pdf",.plot,height = 10,width = 10)


###burden test 同様singleton, doubleton除去して見たら？
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
  dplyr::select(gene_symbol,chr,start,end,ref,alt,focal)
truncating_count_rare_cut_single = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  left_join(truncating_focal_site) %>>%
  filter(!is.na(focal))%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)
truncate_plot_allcantype(truncating_count_rare_cut_single)
}
#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
###### truncating ########
truncating_count_onco_rare = all_maf_for_cumulative %>>%
  filter(MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare.tsv")
ggsave("age_plot/cumulative/onco/truncating_rare.pdf",.plot_all,height = 5,width = 5)
#oncogeneのtruncateあるなしのt_testでは？？
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/onco/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.4394
t.test(truncating_count_onco_rare[truncating_count_onco_raer$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(0){
truncating_count_onco = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

#truncateの数だけでは？？
.plot_all = truncating_count_onco %>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  truncate_plot_allcantype()
.plot_all
ggsave("age_plot/cumulative/onco/truncating.pdf",.plot_all,height = 5,width = 5)
#oncogeneのtruncateあるなしのt_testでは？？
.plot_by = truncating_count_onco %>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  truncate_plot_bycantype()
.plot_by
ggsave("age_plot/cumulative/onco/truncating_by_cancertype.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.4394
t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
       truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/alloncogene.pdf",.plot,width = 10,height = 10)
}
