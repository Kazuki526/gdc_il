######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.R,quality_filter.Rを実行してからこのスクリプトを開始する！！！


#####################################################################################################################
#指数表記
exponent_notation = function(.num){
  .log=trunc(log10(.num))
  .log=ifelse(.log>0,.log+1,.log-1)
  paste0(.num*(10^-.log)," %*% 10^",.log)
}

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
             missense_count_order = missense_num,gender="male",
             age=-50, cancer_type="BRCA", patient_id="no_patient", truncating_count_n=0)
    return(rbind(.missense_count,.outcome))
  }else{return(.missense_count)}
}

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

all_maf_for_cumulative = maf_trim_for_cumulative()

###################################################################################################
##################################################################################################
cumulative_plot = function(.maf=all_maf_for_cumulative,.MAF_start = 0,.mutype="missense",
                           .MAF_end, .path, .role="TSG", .facet_by_cancer_type = F,
                           .by_gene = F, .more_5par = F,.race="all"){
  .patient_list=patient_list
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
  
  #患者ごとのtruncating な遺伝子の数
  .truncating_count = .maf %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter((mutype=="truncating"|mutype=="splice"),role==.role) %>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    #0個の患者も入れる
    {left_join(.patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  if(.mutype=="silent"){
    .maf = .maf %>>%
      group_by(patient_id,gene_symbol,Protein_position) %>>%
      mutate(same_codon=n()) %>>%
      filter(same_codon==1) %>>%
      ungroup()
  }
  ##missense の数
  missense_count = .maf %>>%
    filter(MAF>=.MAF_start/100, MAF<=.MAF_end/100,MAC!=0,mutype==.mutype) %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role) %>>%
    #missenseの数からmisssenseを持つ遺伝子の数に変えるか
    {if(.by_gene){.%>>%group_by(cancer_type,patient_id,gene_symbol) %>>%summarise(MAC=1)}else{.}} %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #MAF>5%だと数が多いから3で割る
    {if(.more_5par){.%>>%mutate(missense_num = missense_num %/% 3)}else{.}} %>>%
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num)) %>>%
    {left_join(.patient_list,.)} %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)　%>>%
    left_join(.truncating_count %>>%dplyr::select(-age)) %>>%
    filter(truncating_count_n==0)
  #相関直線を
  regression=0
  if(.facet_by_cancer_type){  #cancer_typeごとに
    lm_p=function(.data){
      lm=lm(age ~ missense_num, data=.data)
      as.data.frame(as.list(coef(lm))) %>>%
        mutate(p_value_=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                                summary(lm)$fstatistic["dendf"])))%>>%
        mutate(p_value=as.character(ifelse(p_value_<0.001,exponent_notation(signif(p_value_,3)),
                                           signif(p_value_,3))))
    }
    regression = missense_count %>>%
      tidyr::nest(-cancer_type) %>>%
      mutate(data = purrr::map(data,~lm_p(.))) %>>%
      unnest()
  }else{    #全cancer_type
    lm=lm(age ~ missense_num, data=missense_count)
    regression=as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value_=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"])))%>>%
      mutate(p_value=as.character(ifelse(p_value_<0.001,exponent_notation(signif(p_value_,3)),
                                         signif(p_value_,3))))
  }
  write_df(regression,paste0("age_plot/cumulative/",.path,"/lm",.MAF_start,"~",.MAF_end,"regression.tsv"))
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
    ylim(0,90)+
    ggtitle(paste0("MAF = ",.MAF_start," ~ ",.MAF_end," %"))+
    theme(title = element_text(size = 50), axis.title = element_text(size = 50),
          panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray"),
          panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.ticks.y = element_blank())
  if(.facet_by_cancer_type){ #cancer_typeごとに
    .plot=.plot+
      facet_wrap( ~ cancer_type, ncol = 4)+
      geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                  filter(missense_num > 0),aes(x=.coef_posi,y=10,label =coefi),colour = "blue",size=6,hjust=1)+
      geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                  filter(missense_num <= 0),aes(x=.coef_posi,y=10,label =coefi),colour = "red" ,size=6,hjust=1)+
      geom_text(data =missense_count %>>%count(cancer_type,missense_count_n),
                aes(x=missense_count_n,y=5,label=n),size=3,position="stack")+
      geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
                size=7,position = "stack", parse = T)+
      theme(axis.text = element_text(size=16), strip.text = element_text(size=20),
            strip.background = element_rect(fill="transparent", colour = "black"))
  }else{     #全cancer_type
    .plot=.plot+
      geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                  filter(missense_num > 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="blue",size=15,hjust=1)+
      geom_text(data = regression %>>%mutate(coefi=paste0("coef=",round(missense_num*100)/100))%>>%
                  filter(missense_num <= 0),aes(x=.coef_posi,y=10,label =paste0(coefi,"")),colour="red" ,size=15,hjust=1)+
      geom_text(data =missense_count %>>%count(missense_count_n),
                aes(x=missense_count_n,y=5,label=n),size=8,position="stack")+
      geom_text(data = regression,aes(x=.p_posi,y=0,label=paste0("p~ ~value==",p_value)),
                size=20,position = "stack", parse = T)+
      theme(axis.text = element_text(size=40))
  }
  if(.more_5par){  #### MAF>5%の時
    .plot=.plot+scale_x_discrete(labels = legendx3)+
      xlab("number of gene")
  }else{
    .plot = .plot+xlab(paste0("number of ", .mutype, " mutation"))
  }
  .plot
  ggsave(paste0("age_plot/cumulative/",.path,"/plot",.MAF_start,"~",.MAF_end,".pdf"),.plot,height = 16,width = 15)
}
#################################################################################################################
#0.01%ごとにregressionしてみる
make_regression_tabel = function(.maf=all_maf_for_cumulative,.vcf=vcf_exac,.race="all",.role = "TSG",
                                 .fdr=0.01,.mutype="missense",.max_maf=10,.maf_filter=F,
                                 .database="all",.duplicate=T,.somatic=T,.varscan=T){
  regression_out = function(.class,.maf,.role,.patient_list,.truncate){
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .class)%>>%
      left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
                  dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
      filter(role==.role) %>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>%
      {left_join(.patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num)) %>>%
      left_join(.truncate,by = c("cancer_type","patient_id","age")) %>>%
      filter(truncating_count_n==0)
    #相関直線を
    lm=lm(age/365.25 ~ missense_num, data=missense_count)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"]))
  }
  if((.race=="all")&(substitute(.vcf)=="vcf_exac")){
    .vcf = .vcf %>>%
      mutate(AF=AC_Adj/AN_Adj)%>>%
      dplyr::select(chr,start,ref,alt,AF)
  }
  .patient_list=patient_list
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
  if(.maf_filter){.maf=maf_trim_for_cumulative(.vcf=.vcf,.race=.race,.fdr=.fdr,.database=.database,
                               .duplicate=.duplicate,.somatic=.somatic,.varscan=.varscan)
  }
  .truncate=.maf %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter((mutype=="truncating"|mutype=="splice"),role==.role) %>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(.patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  .maf = .maf %>>%filter(mutype==.mutype)
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.role,.patient_list,.truncate)))%>>%
    unnest()
}
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
##############################全cancer_typeまとめて################################################
#ある人とない人で検定してみる
rbind(d_maf,e_maf) %>>% #0.05%以下
  filter(mutype=="missense") %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
  filter(role=="TSG") %>>%
  group_by(patient_id) %>>%
  summarise(focal = "have") %>>%
  {left_join(patient_list%>>%
               left_join(truncating_count)%>>%filter(truncating_count_n==0),.)} %>>%
  mutate(focal = ifelse(is.na(focal),"zero",focal)) %>>%
  {t.test(.[.$focal=="have",]$age/365.25,.[.$focal=="zero",]$age/365.25,alternative = "less")}

###################################################################################################
cumulative_plot(.MAF_end = 0.5,.path = "all_cancer_type")
cumulative_plot(.MAF_end = 0.23,.path = "all_cancer_type")
cumulative_plot(.MAF_end = 0.05,.path = "all_cancer_type")
cumulative_plot(.MAF_start = 0.23,.MAF_end = 0.5,.path = "all_cancer_type")
cumulative_plot(.MAF_start = 0.5,.MAF_end = 1,.path = "all_cancer_type")

#cumulative_plot(.MAF_start = 5,.MAF_end = 50,.path = "all_cancer_type", .more_5par = T)
#cumulative_plot(.MAF_start = 1,.MAF_end = 5,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.5,.MAF_end = 1,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.1,.MAF_end = 0.5,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.05,.MAF_end = 0.1,.path = "all_cancer_type")

####################################################################################################
########################################## cancer_type ごとに ######################################
cumulative_plot(.MAF_end = 0.5,.path = "by_cancer_type",.facet_by_cancer_type = T)
cumulative_plot(.MAF_end = 0.23,.path = "by_cancer_type",.facet_by_cancer_type = T)
cumulative_plot(.MAF_end = 0.05,.path = "by_cancer_type",.facet_by_cancer_type = T)

cumulative_plot(.MAF_start = ,.MAF_end = ,.path = "by_cancer_type",.facet_by_cancer_type = T)
cumulative_plot(.MAF_start = ,.MAF_end = ,.path = "by_cancer_type",.facet_by_cancer_type = T)
cumulative_plot(.MAF_start = ,.MAF_end = ,.path = "by_cancer_type",.facet_by_cancer_type = T)
cumulative_plot(.MAF_start = ,.MAF_end = ,.path = "by_cancer_type",.facet_by_cancer_type = T)

#########################################################################
# mutationのあるgene数でcount
cumulative_plot(.MAF_end = 0.23,.path = "gene",.by_gene = T)
cumulative_plot(.MAF_end = 0.5,.path = "gene",.by_gene = T)
#silent
cumulative_plot(.MAF_end = 0.5,.mutype="silent",.path = "silent")
cumulative_plot(.MAF_end = 0.05,.mutype="silent",.path = "silent")
cumulative_plot(.MAF_end = 0.01,.mutype="silent",.path = "silent")
###################################################################################################
###################################################################################################
regression_tbl_plot = function(.reg_tbl,.maf_max = 1,
                               .bl_ln = NULL,.red_ln = NULL,.expand = 0.03){
  .plot=.reg_tbl　%>>%
    mutate(MAF=MAF*100)%>>%
    filter(MAF<.maf_max)%>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    geom_vline(xintercept = 0,size =1)+
    geom_vline(xintercept = 1, size=0.1)+
    labs(x="MAF (%)",y="regression coefficient")+
    scale_x_continuous(limits = c(0,NA), expand = c(0,.expand))+
    theme_bw()+
    theme(axis.text = element_text(size = 20),axis.title = element_text(size = 30))
  if(.maf_max==10){.plot = .plot + geom_hline(yintercept = 0,size =1)}
  if(!is.null(.bl_ln)) {.plot = .plot + geom_vline(xintercept = .bl_ln,colour="blue")}
  if(!is.null(.red_ln)){.plot = .plot + geom_vline(xintercept = .red_ln,colour="red")}
  plot(.plot)
  return(.plot)
}

regression_table = make_regression_tabel()
.plot=regression_table　%>>%
  regression_tbl_plot(.bl_ln = 0.23,.red_ln = 0.5)
ggsave("age_plot/cumulative/regression_plot-1.pdf",.plot,height = 8,width = 20)
.plot=regression_table　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.23,.red_ln = 0.5,.expand = 0.15)
ggsave("age_plot/cumulative/regression_plot-10.pdf",.plot,height = 8,width = 20)

############suppliment用##############
#silent
regression_table_silent = make_regression_tabel(.mutype = "silent")
.plot= regression_table_silent %>>%
  regression_tbl_plot(.maf_max = 1, .bl_ln = NULL, .red_ln = NULL)
ggsave("age_plot/cumulative/silent/regression_plot-1.pdf",.plot,height = 8,width = 20)
.plot= regression_table_silent %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
ggsave("age_plot/cumulative/silent/regression_plot-10.pdf",.plot,height = 8,width = 20)


#人種白人だけでやってみると？
reg_tbl=make_regression_tabel(.race = "white",.maf_filter = T)
reg_tbl %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
reg_tbl_silent=make_regression_tabel(.mutype = "silent",.race = "white")
reg_tbl_silent%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)

##############################################################################################
#どの遺伝子が1番効いてる？
get_lm_coef = function(.tbl){
  if(first(.tbl$focal =="yes")){
    .patient_count = left_join(patient_list,.tbl, by = "patient_id") %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num),
             age=signif(age/365.25,4))
    model=lm(age ~ missense_num, data = .patient_count)
    as.data.frame(as.list(coef(model))) %>>%
      mutate(lm_p_value = 1 - pf(summary(model)$fstatistic["value"],summary(model)$fstatistic["numdf"],
                                 summary(model)$fstatistic["dendf"]),
             t_test_p_value = t.test(.patient_count[.patient_count$missense_num > 0,]$age/365.25,
                                     .patient_count[.patient_count$missense_num ==0,]$age/365.25,
                                     alternative= "less")$p.value,
             all_missense_num = sum(.patient_count$missense_num))
  }else{
    data.frame(X.Intercept.=NA,missense_num=NA,lm_p_value=NA,t_test_p_value=NA,
               all_missense_num=sum(.tbl$missense_num,na.rm = T))
  }
}
driver_genes %>>%dplyr::select(gene,role) %>>%dplyr::rename(gene_symbol=gene) %>>%
  filter(role=="TSG",gene_symbol!="KMT2C") %>>%
  left_join(all_maf_for_cumlative %>>%filter(MAF<0.005)%>>%
              group_by(gene_symbol,patient_id) %>>% summarise(missense_num=sum(MAC))) %>>%
  mutate(focal = ifelse(is.na(patient_id),"no","yes")) %>>%
  nest(-gene_symbol,-role) %>>%
  mutate(regression =purrr::map(data, ~get_lm_coef(.))) %>>%
  dplyr::select(-role,-data) %>>%unnest() %>>%
  arrange(missense_num) %>>%
  write_df("~/Dropbox/install/tvz/TSG_by_gene_regression_coef.tsv")


#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
####### missense ########
cumulative_plot(.MAF_end = 0.06,.path = "onco",.role = "oncogene")
cumulative_plot(.MAF_end = 3,.path = "onco",.role = "oncogene")

cumulative_plot(.MAF_end = 0.02,.path = "onco/silent",.role = "oncogene",.mutype = "silent")
####### missense by cancer type ########
cumulative_plot(.MAF_end = 0.06,.path = "onco/by_cancer_type",.role = "oncogene",.facet_by_cancer_type = T)
cumulative_plot(.MAF_end = 3,.path = "onco/by_cancer_type",.role = "oncogene",.facet_by_cancer_type = T)

################################################################################################################
regression_table_onco = make_regression_tabel(.role = "oncogene")
.plot=regression_table_onco　%>>%
  regression_tbl_plot(.bl_ln = 0.06)
ggsave("age_plot/cumulative/onco/regression_plot-1.pdf",.plot,height = 8,width = 20)
.plot=regression_table_onco　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.06,.red_ln = 3,.expand = 0.15)
ggsave("age_plot/cumulative/onco/regression_plot-10.pdf",.plot,height = 8,width = 20)

############suppliment用##############
#silent
regression_table_silent_onco = make_regression_tabel(.mutype = "silent",.role = "oncogene")
.plot= regression_table_silent_onco %>>%
  regression_tbl_plot(.maf_max = 1, .bl_ln = NULL, .red_ln = NULL)
ggsave("age_plot/cumulative/onco/silent/regression_plot-1.pdf",.plot,height = 8,width = 20)
.plot= regression_table_silent_onco %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
ggsave("age_plot/cumulative/onco/silent/regression_plot-10.pdf",.plot,height = 8,width = 20)



##################################################################################################################
library(perm)
#MAF>5%で有意に発症年齢を下げているという結果が出たのでサイトごとに見てみる
onco_major=all_maf_for_cumlative%>>%
  filter(MAF>=0.05)%>>%
  left_join(truncating_count_onco) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
  filter(role=="oncogene",truncating_count_n==0) %>>%
  mutate(genotype=ifelse(n_allele2==alt,ifelse(n_allele1==alt,"altalt","refalt"),"refref"),
         genotype_order=ifelse(n_allele2==alt,ifelse(n_allele1==alt,3,2),1),
         site=paste(gene_symbol,chr,start,paste0(ref,">",alt),sep=":"),
         age=signif(age/365.25,4))
#各siteでpermutation testする
site_permP = onco_major %>>%
  nest(-site) %>>%
  mutate(p_value=purrr::map(data,~permKS(.$age,.$genotype,,method = "exact.mc",
                                         control = permControl(nmc=10^4))$p.value)) %>>%
  dplyr::select(-data) %>>% unnest()

.plot = onco_major %>>%
  ggplot(aes(x=reorder(genotype,genotype_order), y=age))+
  geom_violin()+
  geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)+
  facet_wrap( ~ site)+
  ylim(0,95)+
  geom_text(data = onco_major%>>%count(site,genotype),aes(x=genotype,y=0,label=n),size=3)+
  geom_text(data = onco_major%>>%group_by(site,genotype)%>>%summarise(mean_age=signif(mean(age),4)),
            aes(x=genotype,y=10,label=mean_age),colour="red",size=3)+
  geom_text(data = site_permP,aes(x="refalt",y=95,label=ifelse(p_value<0.01,"***"," ")),size=8, colour="red")+
  xlab("geonotype")+
  theme(strip.text = element_text(size=10),panel.background = element_rect(fill="transparent",colour="black"),
        panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
        axis.text = element_text(size=16),panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_line(colour = "gray"),
        axis.ticks.y = element_blank(),axis.title = element_text(size=35),
        strip.background = element_rect(fill="transparent", colour = "black"))
.plot
ggsave("age_plot/onco_5-50_bysite.pdf",.plot,height = 12,width = 10)

################# gene　ごとに？############
driver_genes %>>%dplyr::select(gene,role) %>>%dplyr::rename(gene_symbol=gene) %>>%
  filter(role=="oncogene") %>>%
  left_join(all_maf_for_cumlative %>>%filter(MAF<0.0006)%>>%
              group_by(gene_symbol,patient_id) %>>% summarise(missense_num=sum(MAC))) %>>%
  group_by(gene_symbol)%>>%mutate(focal=ifelse(length(patient_id)>1,"yes","no")) %>>%ungroup()%>>%
  nest(-gene_symbol,-role) %>>%
  mutate(regression = purrr::map(data, ~get_lm_coef(.))) %>>%
  dplyr::select(-role,-data) %>>%unnest() %>>%
  arrange(missense_num) %>>%
  write_df("~/Dropbox/install/tvz/oncogene_by_gene_regression_coef.tsv")
