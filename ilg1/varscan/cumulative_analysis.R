######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.R, quality_filter.Rを実行してからこのスクリプトを開始する！！！


######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
ref_minor = mid_af_coverage %>>%
  left_join(left_join(tally_norm_maf,vcf_exac%>>%
                        mutate(AF=AC_Adj/AN_Adj)%>>%
                        dplyr::select(chr,start,ref,alt,AF))%>>%
              dplyr::select(-ac_cancer,-hom_cancer)) %>>%
  filter(AF > 0.5) %>>% mutate(MAF=1-AF) %>>%
  left_join(norm_maf_all) %>>%
  mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
         n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
         soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
         LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
  dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
                mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,age) %>>%
  quality_filter(.data_type = "maf",.varscan = T) %>>%
  mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
  filter(MAC > 0) %>>%
  dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age)

all_maf_for_cumulative = norm_maf_all %>>%
  quality_filter(.data_type = "maf",.varscan = T) %>>%
  left_join(vcf_exac%>>%mutate(AF=AC_Adj/AN_Adj)%>>%dplyr::select(chr,start,ref,alt,AF)) %>>%
  mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
  filter(AF < 0.5) %>>% mutate(MAF = AF) %>>%
  mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
  dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age) %>>%
  rbind(ref_minor) %>>%
  filter(!is.na(age))



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
truncate_plot_allcantype= function(.tbl){
  .max_count=max(.tbl$truncating_count_n)
  .p_posi=.max_count/2 + 1
  .coef_posi=.max_count+1
  lm=lm(age ~ truncating_count_n, data=.tbl)
  regression=as.data.frame(as.list(coef(lm))) %>>%
    mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                           summary(lm)$fstatistic["dendf"])))%>>%
    mutate(p_value=as.character(ifelse(p_value<0.001,exponent_notation(signif(p_value,3)),
                                       signif(p_value,3))))
  .tbl %>>%truncating_filling()%>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    ylim(0,90)+
    geom_text(data =.tbl %>>%count(truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=5,label=n),size=3,position="stack")+
    geom_abline(data = regression %>>%filter(truncating_count_n > 0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),colour = "blue")+
    geom_abline(data = regression %>>%filter(truncating_count_n <= 0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),colour = "red")+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(truncating_count_n*100)/100))%>>%
              filter(truncating_count_n > 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),
              colour = "blue",size=4)+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(truncating_count_n*100)/100))%>>%
              filter(truncating_count_n <= 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),
              colour = "red",size=4)+
    geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
              size=8,position = "stack", parse = T)+
    xlab("number of truncated gene")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=20), axis.text = element_text(size=20))
}
truncate_plot_bycantype = function(.tbl){
  .max_count=max(.tbl$truncating_count_n)
  .p_posi=.max_count/2 + 1
  .coef_posi=.max_count+1
  lm_p=function(.data){
    lm=lm(age ~ truncating_count_n, data=.data)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                             summary(lm)$fstatistic["dendf"])))%>>%
      mutate(p_value=as.character(ifelse(p_value<0.001,exponent_notation(signif(p_value,3)),
                                         signif(p_value,3))))
  }
  regression = .tbl %>>%
    tidyr::nest(-cancer_type) %>>%
    mutate(data = purrr::map(data,~lm_p(.))) %>>%
    unnest()
  .tbl %>>%truncating_filling()%>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    ylim(0,90)+
    facet_wrap( ~ cancer_type, ncol = 4)+
    geom_text(data =.tbl %>>%count(cancer_type,truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=5,label=n),size=3,position="stack")+
    geom_abline(data = regression %>>%filter(truncating_count_n > 0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),colour = "blue")+
    geom_abline(data = regression %>>%filter(truncating_count_n <= 0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),colour = "red")+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(truncating_count_n*100)/100))%>>%
                filter(truncating_count_n > 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),
              colour = "blue",size=4,hjust=1)+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(truncating_count_n*100)/100))%>>%
                filter(truncating_count_n <= 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),
              colour = "red",size=4,hjust=1)+
    geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
              size=5,position = "stack", parse = T)+
    xlab("number of truncated gene")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=20), axis.text = element_text(size=20),
          strip.background = element_rect(fill="transparent", colour = "black"))
 }
#患者ごとのtruncating な遺伝子の数
truncating_count = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter((mutype=="truncating"|mutype=="splice"),role=="TSG") %>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))


.plot = truncating_count%>>%
  mutate(age=round(age/365.25*100)/100)%>>%
  truncate_plot_allcantype()
.plot
ggsave("age_plot/cumulative/truncating.pdf",.plot,height = 5,width = 5)
.plot = truncating_count%>>%
  mutate(age=round(age/365.25*100)/100)%>>%
  truncate_plot_bycantype()
plot(.plot)
ggsave("age_plot/cumulative/truncating_by_cancerype.pdf",.plot,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p_value=0.02518
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

##### BRCA1,2を別にしてみたら？####
brca_truncate_num = patient_list %>>%
  left_join(all_maf_for_cumulative %>>%
              filter(!(gene_symbol=="BRCA1"|gene_symbol=="BRCA2")) %>>%  ####ここを変えてBRCA1,2のみか以外か変える
              left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
              filter((mutype=="truncating"|mutype=="splice"),role=="TSG") %>>%
              count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
              group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot = brca_truncate_num %>>%truncate_plot_allcantype()
.plot
ggsave("age_plot/cumulative/brca_truncating.pdf",.plot,height = 5,width = 5)
.plot = brca_truncate_num %>>%truncate_plot_bycantype()
.plot
ggsave("age_plot/cumulative/brca_truncating_by_cancerype.pdf",.plot,height = 10,width = 10)

##filterをいじってBRCA1,2以外では？
.plot = brca_truncate_num %>>%truncate_plot_allcantype()
.plot
ggsave("age_plot/cumulative/brca_non_truncating.pdf",.plot,height = 5,width = 5)
.plot = brca_truncate_num %>>%truncate_plot_bycantype()
.plot
ggsave("age_plot/cumulative/brca_non_truncating_by_cancerype.pdf",.plot,height = 10,width = 10)


############### MAF<0.05%のtruncating mutationのみでやってみたら？
truncating_count_rare = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter((mutype=="truncating"|mutype=="splice"),role=="TSG") %>>%
  filter(MAF < 0.05)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)


.plot = truncating_count_rare%>>%truncate_plot_allcantype()
.plot
ggsave("age_plot/cumulative/rare_truncating.pdf",.plot,height = 5,width = 5)
.plot = truncating_count_rare%>>%truncate_plot_bycantype()
plot(.plot)
ggsave("age_plot/cumulative/raer_truncating_by_cancerype.pdf",.plot,height = 10,width = 10)
#あるなしのt.test p-value=0.002675
t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")

#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################

###### truncating ########
truncating_count_onco = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter((mutype=="truncating"|mutype=="splice"),role=="oncogene") %>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

#truncateの数だけでは？？
.plot = truncating_count_onco %>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  truncate_plot_allcantype()
.plot
ggsave("age_plot/cumulative/onco/truncating.pdf",.plot,height = 5,width = 5)
#oncogeneのtruncateあるなしのt_testでは？？
.plot = truncating_count_onco %>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  truncate_plot_bycantype()
.plot
ggsave("age_plot/cumulative/onco/truncating_by_cancertype.pdf",.plot,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.4394
t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
       truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
