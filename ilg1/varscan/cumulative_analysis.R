######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.Rを実行してからこのスクリプトを開始する！！！


######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
######## classわけ (MAF=)####
#(~5, 5~1, 1~0.5, 0.5~0.05, 0.05~)
# a   b    c      d        e
classed_site=tally_norm_maf %>>%
  left_join(vcf_exac)%>>%
  left_join(coverage_all ) %>>%
  mutate(AF = ifelse(is.na(ac_exac),0,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))

mid_maf = mid_af_coverage %>>%
  left_join(classed_site) %>>%
  dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac) %>>%
  filter(AF > 0.005) %>>%
  left_join(somatic_recurrent)%>>%
  filter(is.na(n)) %>>%dplyr::select(-n)%>>%
  remove_duplicate %>>%
  mutate(cancer_type = ifelse((cancer_type=="COAD" | cancer_type=="READ"),"CRC",cancer_type))%>>%
  mutate(cancer_type = ifelse((cancer_type=="KIRC" | cancer_type=="KIRP" | cancer_type=="KICH"),"KCC",cancer_type))%>>%
  left_join(norm_maf %>>%select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                                -Protein_position,-strand,-t_genotype,-n_genotype)%>>%
              filter(!(soma_or_germ =="somatic" & LOH=="no"))%>>%
              mutate(alt=n_allele2)) %>>%
  mutate(t_allele1 = ifelse(is.na(t_allele1),ref,t_allele1),
         t_allele2 = ifelse(is.na(t_allele2),ref,t_allele2),
         n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
         n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
         soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
         LOH = ifelse(is.na(LOH),"ref",LOH))

a_maf = mid_maf %>>%filter(MAF>=0.05, mutype!="silent")
b_maf = mid_maf %>>%filter(MAF<0.05,MAF>=0.01, mutype!="silent")
c_maf = mid_maf %>>%filter(MAF<0.01,MAF>=0.005, mutype!="silent")

ref_minor_ef = mid_maf %>>%
  filter(AF>0.995) %>>%
  filter(ref==n_allele1) %>>%
  mutate(soma_or_germ = ifelse(soma_or_germ == "ref","ref_minor",soma_or_germ))
rm(mid_maf)

d_maf = norm_maf %>>% 
  filter(!(soma_or_germ =="somatic" & LOH=="no"))%>>%
  dplyr::select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                -Protein_position,-strand,-t_genotype,-n_genotype) %>>%
  left_join(classed_site %>>%dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac)) %>>%
  filter(AF<0.005)%>>%
  rbind(ref_minor_ef) %>>%
  filter(MAF<0.005,MAF>=0.0005, mutype!="silent",cancer_type!="KICH")
e_maf = norm_maf %>>%
  filter(!(soma_or_germ =="somatic" & LOH=="no"))%>>%
  dplyr::select(-Transcript_ID,-t_depth,-t_ref_count,-t_alt_count,-n_depth,-n_ref_count,-n_alt_count,
                -Protein_position,-strand,-t_genotype,-n_genotype) %>>%
  left_join(classed_site %>>%dplyr::select(-ac_cancer,-an_cancer,-homo_cancer,-ac_exac,-an_exac)) %>>%
  filter(AF<0.005)%>>%
  rbind(ref_minor_ef) %>>%
  filter(MAF<0.0005, mutype!="silent",cancer_type!="KICH")

#####################################################################################################################
#####################################################################################################################
#violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
missense_filling = function(.missense_count){
  .max_count=max(.missense_count$missense_num)
  .max_count=ifelse(.max_count>10,10,.max_count)
  .outcome = data.frame(missense_num=1:.max_count) %>>%
    left_join(.missense_count%>>%count(missense_num)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$missense_num) != 0){
    .outcome = .outcome %>>%
      mutate(missense_count_n = as.character(missense_num),
           missense_count_order = missense_num,
           age=-50, cancer_type="BRCA", patient_id="no_patient", truncating_count_n=0)
      return(rbind(.missense_count,.outcome))
  }else{return(.missense_count)}
}
truncating_filling = function(.truncating_count){
  .max_count=max(.truncating_count$truncating_count_n)
  .outcome = data.frame(truncating_count_n=1:.max_count) %>>%
    left_join(.truncating_count%>>%count(truncating_count_n)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$truncating_count_n) != 0){
    .outcome = .outcome %>>%
      mutate(age=-50, cancer_type="BRCA", patient_id="no_patient")
    return(rbind(.truncating_count,.outcome))
  }else{return(.truncating_count)}
}
#####################################################################################################################
#指数表記
exponent_notation = function(.num){
  .log=trunc(log10(.num))
  .log=ifelse(.log>0,.log+1,.log-1)
  paste0(.num*(10**-.log),"%*%10^",.log)
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
truncating_count = norm_maf %>>%
  filter(!(soma_or_germ =="somatic" & LOH=="no")) %>>%
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
  left_join(norm_maf %>>%
              filter(!(soma_or_germ =="somatic" & LOH=="no")) %>>%
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
truncating_count_rare = norm_maf %>>%
  filter(!(soma_or_germ =="somatic" & LOH=="no")) %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
  filter((mutype=="truncating"|mutype=="splice"),role=="TSG") %>>%
  left_join(vcf_exac %>>% dplyr::rename(n_allele2=alt) %>>%
              mutate(MAF=ifelse(ac_exac/an_exac < 0.5, ac_exac/an_exac*100, (1-ac_exac/an_exac)*100)))%>>%
  mutate(MAF=ifelse(is.na(MAF),0,MAF))%>>%
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
####################################################################################################################
##############################全cancer_typeまとめて################################################
mid_count = function(.tbl){
  .tbl %>>%
    mutate(alt_count = ifelse(n_allele2==ref,0,ifelse(n_allele1==ref,1,2))) %>>%
    mutate(MAC =ifelse(AF !=MAF,2-alt_count,alt_count)) %>>%
    dplyr::select(-alt_count)
}

cumulative_plot_all = function(.tbl,.class,.role="TSG"){
  ##missense の数
  missense_count =.tbl%>>%
    filter(mutype=="missense") %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role) %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #mutate(missense_num=missense_num %/% 3) %>>% ###############5%以上やる時用
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num))
  if(.role=="TSG"){
    patient_trunc_num=left_join(patient_list,truncating_count)
  }else{
    patient_trunc_num=left_join(patient_list,truncating_count_onco)
  }
  #上のにmissense_count_n=0のpatient情報を追加
  missense_count=patient_trunc_num %>>%
    filter(truncating_count_n==0)%>>%
    left_join(missense_count) %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
  #相関直線を
  lm=lm(age ~ missense_num, data=missense_count)
  regression=as.data.frame(as.list(coef(lm))) %>>%
    mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                           summary(lm)$fstatistic["dendf"])))%>>%
    mutate(p_value=as.character(ifelse(p_value<0.001,exponent_notation(signif(p_value,3)),
                                       signif(p_value,3))))
  if(.role=="TSG"){
    write_df(regression,paste0("age_plot/cumulative/all_cancer_type/lm",.class,"regression.tsv"))
  }else{
    write_df(regression,paste0("age_plot/cumulative/onco/lm",.class,"regression.tsv"))
  }
  legendx3=function(.legend){
    .legend = as.numeric(ifelse(.legend=="10-",10,.legend))
    .legend = .legend*3
    as.character(ifelse(.legend==30,"30-",.legend))
  }
  .max_count=max(missense_count$missense_count_order)
  .p_posi=.max_count/2 +1
  .coef_posi=.max_count +1
  ##バイオリンプロットで見やすく
  .plot = missense_count %>>%missense_filling()%>>%
    ggplot(aes(x=reorder(missense_count_n,missense_count_order), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    geom_abline(data = regression %>>%filter(missense_num > 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "blue")+
    geom_abline(data = regression %>>%filter(missense_num <= 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "red")+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num > 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="blue",size=15,hjust=1)+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num <= 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="red" ,size=15,hjust=1)+
    ylim(0,90)+
    geom_text(data =missense_count %>>%count(missense_count_n),
              aes(x=missense_count_n,y=5,label=n),size=8,position="stack")+
    geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
              size=20,position = "stack", parse = T)+
    xlab(paste0("MAF = ",.class,"%"))+
    #scale_x_discrete(labels = legendx3)+ ####################5%以上のときのみ
    theme(axis.text = element_text(size=50),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.ticks.y = element_blank(),panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray"),
          axis.title.x = element_text(size = 50),axis.title.y = element_blank())
  if(.role=="TSG"){
    ggsave(paste0("age_plot/cumulative/all_cancer_type/plot",.class,".pdf"),.plot,height = 16,width = 15)
  }else{
    ggsave(paste0("age_plot/cumulative/onco/plot",.class,".pdf"),.plot,height = 16,width = 15)
  }
}


#a_maf %>>%mid_count%>>%cumulative_plot_all("5~50") #(この場合はcount数の表示をいじる必要がある)
b_maf %>>%mid_count%>>%cumulative_plot_all("1~5")
c_maf %>>%mid_count%>>%cumulative_plot_all("0.5~1")
#d_maf %>>%mutate(MAC=1) %>>%cumulative_plot_all("plot0.05-0.5")
d_maf %>>%filter(MAF>=0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("0.1~0.5")
d_maf %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("0.05~0.1")
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot_all("0~0.05")
#rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%cumulative_plot_all("plot-0.5")
#rbind(d_maf,e_maf) %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("plot-0.1")
#rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%rbind(c_maf%>>%mid_count())%>>%cumulative_plot_all("plot-1")
#rbind(b_maf,c_maf)%>>%mid_count()%>>%
#    rbind(rbind(d_maf,e_maf)%>>%mutate(MAC=1))%>>%cumulative_plot_all("plot-5")

####################################################################################################
########################################## cancer_type ごとに ######################################
cumulative_plot = function(.tbl,.class,.role="TSG"){
  ##missense の数
  missense_count = .tbl %>>%
    filter(mutype=="missense") %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role) %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #mutate(missense_num=missense_num %/% 3) %>>% ################a_mafをやる時用
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num))
  if(.role=="TSG"){
    patient_trunc_num=left_join(patient_list,truncating_count)
  }else{
    patient_trunc_num=left_join(patient_list,truncating_count_onco)
  }
  #上のにmissense_count_n=0のpatient情報を追加
  missense_count=patient_trunc_num %>>%
    filter(truncating_count_n==0)%>>%
    left_join(missense_count) %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
  #相関直線を
  lm_p=function(.data){
    lm=lm(age ~ missense_num, data=.data)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                             summary(lm)$fstatistic["dendf"])))%>>%
      mutate(p_value=as.character(ifelse(p_value<0.001,exponent_notation(signif(p_value,3)),
                                         signif(p_value,3))))
  }
  regression = missense_count %>>%
    tidyr::nest(-cancer_type) %>>%
    mutate(data = purrr::map(data,~lm_p(.))) %>>%
    unnest()
  if(.role=="TSG"){
    write_df(regression,paste0("age_plot/cumulative/by_cancer_type/lm",.class,"regression.tsv"))
  }else{
    write_df(regression,paste0("age_plot/cumulative/onco/by_cancer_type/lm",.class,"regression.tsv"))
  }
  ##バイオリンプロットで見やすく
  legendx3=function(.legend){
    .legend = as.numeric(ifelse(.legend=="10-",10,.legend))
    .legend = .legend*3
    as.character(ifelse(.legend==30,"30-",.legend))
  }
  .max_count=max(missense_count$missense_count_order)
  .p_posi=.max_count/2 +1
  .coef_posi=.max_count +1
  .plot = missense_count %>>%missense_filling()%>>%
    ggplot(aes(x=reorder(missense_count_n,missense_count_order), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    facet_wrap( ~ cancer_type, ncol = 4)+
    geom_abline(data = regression %>>%filter(missense_num > 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "blue")+
    geom_abline(data = regression %>>%filter(missense_num <= 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "red")+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num > 0),aes(x=.coef_posi,y=10,label =coefi),colour = "blue",size=6,hjust=1)+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num <= 0),aes(x=.coef_posi,y=10,label =coefi),colour = "red" ,size=6,hjust=1)+
    ylim(0,90)+
    geom_text(data =missense_count %>>%count(cancer_type,missense_count_n),
              aes(x=missense_count_n,y=5,label=n),size=3,position="stack")+
    geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
              size=7,position = "stack", parse = T)+
    xlab(paste0("number of missense mutation (MAF = ",.class,"%)"))+
    #scale_x_discrete(labels = legendx3)+ ####################5%以上のときのみ
    theme(strip.text = element_text(size=20),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.text = element_text(size=16),panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray"),
          axis.ticks.y = element_blank(),axis.title = element_text(size=35),
          strip.background = element_rect(fill="transparent", colour = "black"))
  .plot
  if(.role=="TSG"){
    ggsave(paste0("age_plot/cumulative/by_cancer_type/plot",.class,".pdf"),.plot,height = 15,width = 15)
  }else{
    ggsave(paste0("age_plot/cumulative/onco/by_cancer_type/plot",.class,".pdf"),.plot,height = 15,width = 15)
  }
}

#a_maf %>>%mid_count%>>%cumulative_plot("a") #(この場合はcount数の表示をいじる必要がある)
#b_maf %>>%mid_count%>>%cumulative_plot("plot1-5")
#c_maf %>>%mid_count%>>%cumulative_plot("plot0.5-1")
#d_maf %>>%filter(MAF>=0.001)%>>%mutate(MAC=1) %>>%cumulative_plot("plot0.1-0.5")
#d_maf %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot("plot0.05-0.1")
#d_maf %>>%mutate(MAC=1) %>>%cumulative_plot("plot0.05-0.5")
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot("0~0.05")
e_maf %>>%mutate(MAC=1) %>>%filter(MAF<0.0002)%>>%cumulative_plot("0~0.02")
#rbind(b_maf,c_maf)%>>%mid_count()%>>%cumulative_plot("bc")
#c_maf %>>%mid_count()%>>%rbind(d_mid %>>%mutate(MAC=1))%>>%cumulative_plot("cd")
rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%cumulative_plot("0~0.5")
rbind(d_maf,e_maf) %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot("0~0.1")
#rbind(b_maf,c_maf)%>>%mid_count()%>>%rbind(d_mid %>>%mutate(MAC=1))%>>%cumulative_plot("bcd")
rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%rbind(c_maf%>>%mid_count())%>>%cumulative_plot("0~1")
#rbind(b_maf,c_maf)%>>%mid_count()%>>%
#  rbind(rbind(d_maf,e_maf)%>>%mutate(MAC=1))%>>%cumulative_plot("bcde")


#########################################################################
# mutationのあるgene数でcount
cumulative_plot_gene = function(.tbl,.class,.role="TSG"){
  ##missense の数
  missense_count = .tbl %>>%
    filter(MAC>0,mutype=="missense")%>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role) %>>%
    group_by(cancer_type,patient_id,gene_symbol) %>>%summarise(MAC=1) %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #mutate(missense_num=missense_num %/% 2) %>>% ###############5%以上やる時用
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num))
  if(.role=="TSG"){
    patient_trunc_num=left_join(patient_list,truncating_count)
  }else{
    patient_trunc_num=left_join(patient_list,truncating_count_onco)
  }
  #上のにmissense_count_n=0のpatient情報を追加
  missense_count=patient_trunc_num %>>%
    filter(truncating_count_n==0)%>>%
    left_join(missense_count) %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
  #相関直線を
  lm=lm(age ~ missense_num, data=missense_count)
  regression=as.data.frame(as.list(coef(lm))) %>>%
    mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                           summary(lm)$fstatistic["dendf"])))%>>%
    mutate(p_value=as.character(ifelse(p_value<0.001,exponent_notation(signif(p_value,3)),
                                       signif(p_value,3))))
  if(.role=="TSG"){
    write_df(regression,paste0("age_plot/cumulative/gene/lm",.class,"regression.tsv"))
  }else{
    write_df(regression,paste0("age_plot/cumulative/onco/gene/lm",.class,"regression.tsv"))
  }
  legendx2=function(.legend){
    .legend = as.numeric(ifelse(.legend=="10-",10,.legend))
    .legend = .legend*2
    as.character(ifelse(.legend==20,"20-",.legend))
  }
  .max_count=max(missense_count$missense_count_order)
  .p_posi=.max_count/2 +1
  .coef_posi=.max_count +1
  ##バイオリンプロットで見やすく
  .plot = missense_count %>>%missense_filling()%>>%
    ggplot(aes(x=reorder(missense_count_n,missense_count_order), y=age))+
    geom_violin()+
    geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
    geom_abline(data = regression %>>%filter(missense_num > 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "blue")+
    geom_abline(data = regression %>>%filter(missense_num <= 0),
                aes(intercept = X.Intercept.,slope = missense_num),colour = "red")+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num > 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="blue",size=15,hjust=1)+
    geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                filter(missense_num <= 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="red" ,size=15,hjust=1)+
    ylim(0,90)+
    geom_text(data =missense_count %>>%count(missense_count_n),
              aes(x=missense_count_n,y=5,label=n),size=8,position="stack")+
    geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
              size=20,position = "stack", parse = T)+
    xlab(paste0("MAF = ",.class,"%"))+
    #scale_x_discrete(labels = legendx2)+ ####################5%以上のときのみ
    theme(axis.text = element_text(size=50),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray"),axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 50))
  #plot(.plot)
  if(.role=="TSG"){
    ggsave(paste0("age_plot/cumulative/gene/plot",.class,".pdf"),.plot,height = 15,width = 15)
  }else{
    ggsave(paste0("age_plot/cumulative/onco/gene/plot",.class,".pdf"),.plot,height = 15,width = 15)
  }
}

#a_maf %>>%mid_count%>>%cumulative_plot_gene("5~50") #(この場合はcount数の表示をいじる必要がある)
b_maf %>>%mid_count%>>%cumulative_plot_gene("1~5")
c_maf %>>%mid_count%>>%cumulative_plot_gene("0.5~1")
d_maf %>>%filter(MAF>=0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_gene("0.1~0.5")
d_maf %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_gene("0.05~0.1")
#d_maf %>>%mutate(MAC=1) %>>%cumulative_plot_gene("plot0.05-0.5")
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot_gene("0~0.05")

#rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%cumulative_plot_gene("plot-0.5")
#rbind(d_maf,e_maf) %>>%filter(MAF<000.1)%>>%mutate(MAC=1) %>>%cumulative_plot_gene("plot-0.1")


#################################################################################################################
all_site =rbind(a_maf,b_maf)%>>%rbind(c_maf)%>>%mid_count()%>>%
  rbind(rbind(d_maf,e_maf)%>>%mutate(MAC=1))%>>%
  filter(MAC!=0)
regression_out = function(.class,.role="TSG"){
  ##missense の数
  missense_count = all_site %>>%
    filter(MAF<.class)%>>%
    filter(mutype=="missense") %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role) %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC))
  if(.role=="TSG"){
    patient_trunc_num=left_join(patient_list,truncating_count)
  }else{
    patient_trunc_num=left_join(patient_list,truncating_count_onco)
  }
  #上のにmissense_count_n=0のpatient情報を追加
  missense_count=patient_trunc_num %>>%
    filter(truncating_count_n==0)%>>%
    left_join(missense_count) %>>%
    mutate(missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
  #相関直線を
  lm=lm(age ~ missense_num, data=missense_count)
  as.data.frame(as.list(coef(lm)))
}

regression_table = data.frame(MAF=1:1000) %>>%
    mutate(MAF = MAF/10000) %>>%#head(5)%>>%
    mutate(regression = purrr::map(MAF,~regression_out(.)))%>>%
    unnest()
.plot=regression_table　%>>%
  mutate(MAF=MAF*100)%>>%
  #filter(MAF<1)%>>%
  ggplot()+
  geom_point(aes(x=MAF,y=missense_num))+
  geom_vline(xintercept = c(0,0.5),colour="red")+
  geom_vline(xintercept = 1, size=0.1)+
  labs(x="MAF (%)",y="regression coefficient")+
  scale_x_continuous(limits = c(0,NA), expand = c(0,0.1))+#MAF<1の時0.01,MAF<10の時0.1に
  theme_bw()+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 30))
.plot
ggsave("age_plot/cumulative/regression_plot-1.pdf",.plot,height = 8,width = 20)
ggsave("age_plot/cumulative/regression_plot-10.pdf",.plot,height = 8,width = 20)


#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
###### truncating ########
truncating_count_onco = norm_maf %>>%
  filter(!(soma_or_germ =="somatic" & LOH=="no")) %>>%
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

####### missense ########
a_maf %>>%mid_count%>>%cumulative_plot_all("5~50",.role = "oncogene") ##oncoはそれほどmutationないので3で悪必要なし
b_maf %>>%mid_count%>>%cumulative_plot_all("1~5",.role="oncogene")
c_maf %>>%mid_count%>>%cumulative_plot_all("0.5~1",.role="oncogene")
#d_maf %>>%mutate(MAC=1) %>>%cumulative_plot_all("plot0.05-0.5",.role="oncogene")
d_maf %>>%filter(MAF>=0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("0.1~0.5",.role="oncogene")
d_maf %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("0.05~0.1",.role="oncogene")
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot_all("0~0.05",.role="oncogene")
#rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%cumulative_plot_all("plot-0.5",.role="oncogene")
#rbind(d_maf,e_maf) %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot_all("plot-0.1",.role="oncogene")
#rbind(d_maf,e_maf) %>>%mutate(MAC=1) %>>%rbind(c_maf%>>%mid_count())%>>%cumulative_plot_all("plot-1",.role="oncogene")
####### missense by cancer type ########
e_maf %>>%mutate(MAC=1) %>>%cumulative_plot("0~0.05",.role="oncogene")
rbind(d_maf,e_maf) %>>%filter(MAF<0.001)%>>%mutate(MAC=1) %>>%cumulative_plot("0~0.1",.role="oncogene")
d_maf %>>%filter(MAF>=0.001)%>>%mutate(MAC=1) %>>%cumulative_plot("0.1~0.5",.role="oncogene")
c_maf %>>%mid_count%>>%cumulative_plot("0.5~1",.role="oncogene")

################################################################################################################
regression_table_onco =　data.frame(MAF=1:1000) %>>%
  mutate(MAF = MAF/10000) %>>%#head(5)%>>%
  mutate(regression = purrr::map(MAF,~regression_out(.,"oncogene")))%>>%
  unnest()
.plot=regression_table_onco　%>>%
  mutate(MAF=MAF*100)%>>%
  #filter(MAF<1)%>>%
  ggplot()+
  geom_point(aes(x=MAF,y=missense_num))+
  geom_vline(xintercept = c(0,0.1),colour="red")+
  geom_vline(xintercept = 1, size=0.1)+
  labs(x="MAF (%)",y="regression coefficient")+
  scale_x_continuous(limits = c(0,NA), expand = c(0,0.1))+
  theme_bw()+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 30))
.plot
ggsave("age_plot/cumulative/onco/regression_plot-1.pdf",.plot,height = 8,width = 20)
ggsave("age_plot/cumulative/onco/regression_plot-10.pdf",.plot,height = 8,width = 20)
