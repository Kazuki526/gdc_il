###################################################################################################
##################################################################################################
cumulative_plot = function(.maf=all_maf_for_cumulative,.MAF_end, .path="",.title=F,
                           .MAF_start = 0,.mutype="missense",.role="TSG",
                           .facet_by_cancer_type = F, .by_gene = F, .more_5par = F,.race="all",
                           .regression_size = 7,.pnum_size = 4,.save = T,
                           .patient_list=patient_list,.height=8,.width=8,
                           .permu=T,.test_tail="one",.permu_file=NA,.all_color="black"){
  #指数表記
  exponent_notation = function(.num){
    .log=trunc(log10(.num))
    .log=ifelse(.log>0,.log+1,.log-1)
    paste0(.num*(10^-.log)," %*% 10^",.log)
  }
  perm_pvalue = function(.tbl,.regression,.tail,.CT=NA){
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
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
  
  #患者ごとのtruncating な遺伝子の数
  .truncating_count = .maf %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role==.role | role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    #0個の患者も入れる
    {left_join(.patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  ##missense の数
  missense_count = .maf %>>%
    filter(MAF>=.MAF_start/100, MAF<=.MAF_end/100,MAC!=0,mutype==.mutype) %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
    filter(role==.role | role=="oncogene/TSG")%>>%
    #missenseの数からmisssenseを持つ遺伝子の数に変えるか
    {if(.by_gene){.%>>%group_by(cancer_type,patient_id,gene_symbol) %>>%summarise(MAC=1)}else{.}} %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #MAF>5%だと数が多いから3で割る
    {if(.more_5par){.%>>%mutate(missense_num = missense_num %/% 3)}else{.}} %>>%
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num)) %>>%
           {left_join(.patient_list,.)} %>>%
    left_join(patient_with_ps) %>>%filter(is.na(significance))%>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)　%>>%
    left_join(.truncating_count %>>%dplyr::select(-age)) %>>%
    filter(truncating_count_n==0)
  #相関直線を
  regression=0
  if(.facet_by_cancer_type){  #cancer_typeごとに
    lm_p=function(.data,.permu){
      lm=lm(age ~ missense_num, data=.data)
      .dtbl = as.data.frame(as.list(coef(lm)))
      if(.permu){return(.dtbl)
      }else{
        return(mutate(.dtbl,p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                                            summary(lm)$fstatistic["dendf"]))))
      }
    }
    snp_permu_byCT = function(.times,.tbl){
      if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
      .tbl%>>%mutate(age=sample(age,length(age)))%>>%
        nest(-cancer_type)%>>%dplyr::rename(data_=data)%>>%
        mutate(data_ = purrr::map(data_,~lm_p(.,.permu = T)))%>>%unnest()
    }
    regression = missense_count %>>%
      tidyr::nest(-cancer_type) %>>%
      mutate(data = purrr::map(data,~lm_p(.,.permu = .permu))) %>>%
      unnest()
    if(.permu){#if do permutation
      if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
        regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
      }else{
        regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
          mutate(data = purrr::map(times,~snp_permu_byCT(.times =.,.tbl = missense_count))) %>>%
          unnest() %>>%
          rename(regression = missense_num)
        write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
      }
      #regression %>>%
      # mutate(p_value = pmap(.,~ perm_pvalue(.CT = cancer_type,.regression = missense_num,
      #                           .tail = "one",.tbl = regression_tbl))) %>>%unnest()
      regression = regression %>>%rowwise() %>>%
        mutate(p_value = perm_pvalue(.CT = cancer_type,.regression = missense_num,
                                     .tail = .test_tail,.tbl = regression_tbl)) %>>%(?.)
    }
  }else{    #全cancer_type
    snp_permu = function(.tbl,.times,.print=T){
      if(.times %% 1000 == 0 & .print) {print(paste0("permutation ",.times," times now"))}
      .tbl_sample=.tbl%>>%mutate(age=sample(age,length(age)))
      as.tibble(as.list(coef(lm(age/365.25 ~ missense_num,data = .tbl_sample))))
    }
    lm=lm(age ~ missense_num, data=missense_count)
    regression=as.data.frame(as.list(coef(lm)))
    if(.permu){
      if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
        regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
      }else{
        regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
          mutate(tbl=purrr::map(times,~snp_permu(.times = .,.tbl = missense_count))) %>>%
          unnest() %>>%
          rename(regression = missense_num)
        write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
      }
      .p=perm_pvalue(regression_tbl,regression$missense_num,.test_tail)
      regression = mutate(regression,p_value = .p,cancer_type = "All Cancer Types")%>>%(?.)
    }else{
      regression %>>%
        mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                               summary(lm)$fstatistic["dendf"]))) %>>%(?.)%>>%
        mutate(cancer_type="All Cancer Types")
    }
  }
  if(.save ){
    write_df(regression,
             paste0("age_plot/cumulative/",.path,"/lm",.MAF_start,"~",.MAF_end,"regression.tsv"))
  }
  legendx3=function(.legend){
    .legend = as.numeric(ifelse(.legend=="10-",10,.legend))
    .legend = .legend*3
    as.character(ifelse(.legend==30,"30-",.legend))
  }
  regression=regression %>>%
    mutate(out_reg=paste0("R=",signif(missense_num,2),", P",
                          ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  #violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
  missense_filling = function(.missense_count){
    .max_count=max(.missense_count$missense_num)
    .max_count=ifelse(.max_count>9,9,.max_count)
    .outcome = data.frame(missense_num=1:.max_count) %>>%
      left_join(.missense_count%>>%count(missense_num)) %>>%
      filter(is.na(n)) %>>%dplyr::select(-n)
    if(length(.outcome$missense_num) != 0){
      .outcome = .outcome %>>%
        mutate(missense_count_n = as.character(missense_num),
               missense_count_order = missense_num,gender="male",
               age=-50, cancer_type="BRCA", patient_id="no_patient", truncating_count_n=0,
               significance=NA)
      return(rbind(.missense_count,.outcome))
    }else{return(.missense_count)}
  }
  .max_count=max(missense_count$missense_count_order)
  .p_posi=.max_count/2 +1
  .coef_posi=.max_count +1
  ##バイオリンプロットで見やすく
  .plot = missense_count %>>%missense_filling()%>>%
    mutate(cancer_type = ifelse(rep(.facet_by_cancer_type,length(cancer_type)), cancer_type,"All Cancer Types"))%>>%
    ggplot(aes(x=reorder(missense_count_n,missense_count_order), y=age))+
    geom_violin(scale = "count",colour=.all_color)+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill=.all_color,colour=.all_color,shape=21,size=2)+
    geom_abline(data = regression %>>%filter(p_value >= 0.05),
                aes(intercept = X.Intercept.,slope = missense_num),linetype="dashed",colour=.all_color)+
    geom_abline(data = regression %>>%filter(p_value < 0.05),
                aes(intercept = X.Intercept.,slope = missense_num),colour=.all_color)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    ggtitle(if(.title){if(.MAF_start==0){paste0("MAF < ",.MAF_end," %")
    }else{paste0("MAF = ",.MAF_start," ~ ",.MAF_end," %")}}else{""})+
    theme(title = element_text(size = 20), axis.title = element_text(size = 15),
          panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray85"),
          panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.ticks.y = element_blank())
  if(.facet_by_cancer_type){ #cancer_typeごとに
    .plot=.plot+
      facet_wrap( ~ cancer_type, ncol = 5)+
      geom_text(data = regression,aes(x=.coef_posi,y=15,label =out_reg),
                size=.regression_size,hjust=1,colour=.all_color)+
      geom_text(data =missense_count %>>%count(cancer_type,missense_count_n),
                aes(x=missense_count_n,y=5,label=n)
                ,size=.pnum_size,position="stack",angle=-45,colour=.all_color)+
      theme(axis.text = element_text(size=10), strip.text = element_text(size=12),
            strip.background = element_rect(fill="transparent", colour = "black"))
  }else{     #全cancer_type
    .plot=.plot+
      facet_wrap( ~ cancer_type)+
      geom_text(data = regression,aes(x=.coef_posi,y=15,label =out_reg),
                size=.regression_size,hjust=1,colour=.all_color)+
      geom_text(data =missense_count %>>%count(missense_count_n),
                aes(x=missense_count_n,y=5,label=n),size=.pnum_size,position="stack",
                colour=.all_color)+
      theme(axis.text = element_text(size=15), strip.text = element_text(size=12),
            strip.background = element_rect(fill="transparent", colour = "black"))
  }
  .ns = ifelse(.mutype == "missense", "Nonsynonymous","Synonymous")
  if(.more_5par){  #### MAF>5%の時
    .plot=.plot+scale_x_discrete(labels = legendx3)
  }
  if(.by_gene){ ### mutation数でなくmutationを持つ遺伝子数の場合
    .plot = .plot+xlab(paste0("number of gene having ",.ns," Variants"))+
      theme(axis.title = element_text(size=20))
  }else{
    .plot = .plot+xlab(paste0("Number of ", .ns, " Variants"))
  }
  .plot = .plot+ylab("Age at Diagnosis")+guides(colour=F)
  plot(.plot)
  if(.save){
    ggsave(paste0("age_plot/cumulative/",.path,"/plot",.MAF_start,"~",.MAF_end,".pdf"),
           .plot,height = .height,width = .width)
  }
  return(.plot)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#0.01%ごとにregressionしてみる
make_regression_tabel = function(.maf=all_maf_for_cumulative,.vcf=vcf_exac,.race="all",.role = "TSG",
                                 .fdr=0.01,.mutype="missense",.max_maf=10,.maf_filter=F,
                                 .database="all",.duplicate=T,.somatic=T,.varscan=T,
                                 .patient_list = patient_list,.remove0=F){
  regression_out = function(.class,.maf,.role,.patient_list){
    if((.class*10000) %% 100 == 0){print(paste0("doing MAF=",.class*100))}
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .class)%>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>% ungroup() %>>%
      {left_join(.patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num))
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
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
    if(.race == "white"){
      .vcf = .vcf %>>%
        mutate(AF = (AC_FIN + AC_NFE)/(AN_FIN + AN_NFE)) %>>%
        dplyr::select(chr,start,ref,alt,AF,AC_FIN,AC_NFE,AN_FIN,AN_NFE)
    }else{
      .vcf = .vcf %>>%
        mutate(AF = AC_AFR/AN_AFR) %>>%
        dplyr::select(chr,start,ref,alt,AF,AC_AFR,AN_AFR)
    }
  }
  if(.maf_filter){.maf=maf_trim_for_cumulative(.vcf=.vcf,.race=.race,.fdr=.fdr,.database=.database,
                                               .duplicate=.duplicate,.somatic=.somatic,.varscan=.varscan)
  }
  .patient_list=.maf %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role==.role | role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(.patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
    left_join(patient_with_ps,by=c("cancer_type","patient_id")) %>>%filter(is.na(significance))%>>%
    dplyr::select(-significance,-truncating_count_n)
  .maf = .maf %>>%filter(mutype==.mutype) %>>%
    filter(if(.remove0){MAF!=0}else{chr==chr}) %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
                dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
    filter(role==.role | role=="oncogene/TSG")
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.role,.patient_list)))%>>%
    unnest()
}
#### regressionの図のplot 3pattern
regression_plot_log = function(.reg_tbl,.max_maf=50,.min=NA,.dred=NA,.blue=NA,.green=NA){
  .max=max(.reg_tbl$missense_num)
  .min=ifelse(is.na(.min),min(.reg_tbl$missense_num),.min)
  if(.max <0){.max=0}else if(.min >0){.min=0}
  .expand = (.max - .min)/10
  .plot = .reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    #geom_vline(xintercept = 0,size =1)+
    #geom_hline(yintercept = 0,size =1)+
    labs(y="Regression Coefficient",x = "MAF (%)")+
    scale_x_log10(limits = c(0.01,.max_maf), expand = c(0,0.03),
                  breaks=c(0.01,0.1,1,10,50))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,.expand))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          panel.border = element_blank(),axis.line = element_line())
  if(!is.na(.dred)){.plot=.plot+geom_vline(xintercept = .dred,colour="darkred")}
  if(!is.na(.blue)){.plot=.plot+geom_vline(xintercept = .blue,colour="blue")}
  if(!is.na(.green)){.plot=.plot+geom_vline(xintercept = .green,colour="green")}
  plot(.plot)
  return(.plot)
}




regression_plot_byside = function(.reg_tbl,.cut_MAF=1,.max_maf=50,.no_axis_title=F,
                                  .min=NA,.red=0.05,.blue=NA,.green=NA,
                                  .x_axis="yes",.y_axis="yes",.axis_size=12){
  .max=max(.reg_tbl$missense_num)
  .min=ifelse(is.na(.min),min(.reg_tbl$missense_num),.min)
  if(.max <0){.max=0}else if(.min >0){.min=0}
  .expand = (.max - .min)/20
  .plot1 = .reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    filter(MAF<=.cut_MAF) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    #geom_vline(xintercept = 0,size =1)+
    #geom_hline(yintercept = 0,size =1)+
    labs(y="Regression Coefficient")+
    scale_x_continuous(limits = c(0,.cut_MAF), expand = c(0,0.01))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,.expand))+
    theme_bw()+
    theme(axis.text = element_text(size = .axis_size),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_blank(),
          panel.border = element_blank(),axis.line = element_line())
  .plot2 = .reg_tbl%>>%
    mutate(MAF=100*MAF) %>>%
    filter(MAF>=.cut_MAF,MAF<=.max_maf) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    #geom_hline(yintercept = 0,size =1)+
    scale_x_continuous(limits = c(.cut_MAF,.max_maf), expand = c(0,0.1))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,.expand))+
    theme_bw()+
    theme(axis.text.x = element_text(size = .axis_size),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title = element_blank(),panel.border = element_blank(),
          axis.line.x = element_line(colour = "black"))
  if(!is.na(.red) & .red<=.cut_MAF){.plot1=.plot1+geom_vline(xintercept = .red,colour="red")
  }else if(!is.na(.red) & .red<.max_maf){.plot2=.plot2+geom_vline(xintercept = .red,colour="red")}
  if(!is.na(.blue) & .blue<=.cut_MAF){.plot1=.plot1+geom_vline(xintercept = .blue,colour="blue")
  }else if(!is.na(.blue) & .blue<.max_maf){.plot2=.plot2+geom_vline(xintercept = .blue,colour="blue")}
  if(!is.na(.green) & .green<=.cut_MAF){.plot1=.plot1+geom_vline(xintercept = .green,colour="green")
  }else if(!is.na(.green) & .green<.max_maf){.plot2=.plot2+geom_vline(xintercept = .green,colour="green")}
  if(.x_axis=="no"){
    .plot1=.plot1+theme(axis.text.x = element_blank())
    .plot2=.plot2+theme(axis.text.x = element_blank(),axis.line.x = element_line(colour = "black"))}
  if(.y_axis=="no"){.plot1=.plot1+theme(axis.text.y = element_blank())}
  if(.no_axis_title){
    .plot = cowplot::ggdraw()+
      cowplot::draw_plot(.plot2,x=0.64,y=0.05,width = 0.34,height = 0.95)+
      cowplot::draw_plot(.plot1+theme(axis.title = element_blank()),
                                      x=0,y=0.05,width = 0.64,height = 0.95)
  } else {
    .plot = cowplot::ggdraw()+
      cowplot::draw_plot(.plot2,x=0.64,y=0.05,width = 0.34,height = 0.95)+
      cowplot::draw_plot(.plot1,x=0,y=0.05,width = 0.65,height = 0.95)+
      cowplot::draw_text("MAF (%)",x=0.5,y=0.01,vjust = 0,size = 15)
  }
  plot(.plot)
  return(.plot)
}



regression_plot_in = function(.reg_tbl,.in){
  .max=max(.reg_tbl$missense_num)
  .min=min(.reg_tbl$missense_num)
  if(.max <0){.max=0}else if(.min >0){.min=0}
  .expand = (.max - .min)/20
  .max=.max+.expand
  .min=.min-.expand
  .plot1 = .reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    filter(MAF<=1) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    geom_vline(xintercept = 0,size =1)+
    geom_vline(xintercept = 0.05,colour="blue")+
    labs(y="Regression Coefficient",x="MAF (%)")+
    scale_x_continuous(limits = c(0,1), expand = c(0,0.01))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15))
  .plot10=.reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    geom_vline(xintercept = 0,size =1)+
    geom_hline(yintercept = 0,size =1)+
    geom_vline(xintercept = 0.05,colour="blue")+
    labs(y="Regression Coefficient",x="MAF (%)")+
    annotate("rect",xmin=0,xmax=1,ymin=.min,ymax=.max,alpha=0.2)+
    scale_x_continuous(limits = c(0,10), expand = c(0,0.15))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,0))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15))
  if(.in==1){
    .plot=cowplot::ggdraw()+
      cowplot::draw_plot(.plot10,x=0,y=0,width = 1,height = 1)+
      cowplot::draw_plot(.plot1+theme(axis.title = element_blank()),
                         x=0.45,y=0.15,width = 0.5,height = 0.5)
  }else if(.in==10){
    .plot=cowplot::ggdraw()+
      cowplot::draw_plot(.plot1,x=0,y=0,width = 1,height = 1)+
      cowplot::draw_plot(.plot10+theme(axis.title = element_blank()),
                         x=0.45,y=0.15,width = 0.5,height = 0.5)
  }
  plot(.plot)
  return(.plot)
}












#########################################################################
#昔の
if(0){
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
    theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15))
  if(.maf_max==10){.plot = .plot + geom_hline(yintercept = 0,size =1)}
  if(!is.null(.bl_ln)) {.plot = .plot + geom_vline(xintercept = .bl_ln,colour="blue")}
  if(!is.null(.red_ln)){.plot = .plot + geom_vline(xintercept = .red_ln,colour="red")}
  plot(.plot)
  return(.plot)
}
regression_plot_in = function(.reg_tbl,.in){
  .max=max(.reg_tbl$missense_num)
  .min=min(.reg_tbl$missense_num)
  if(.max <0){.max=0}else if(.min >0){.min=0}
  .expand = (.max - .min)/20
  .max=.max+.expand
  .min=.min-.expand
  .plot1 = .reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    filter(MAF<=1) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    geom_vline(xintercept = 0,size =1)+
    geom_vline(xintercept = 0.05,colour="blue")+
    labs(y="Regression Coefficient",x="MAF (%)")+
    scale_x_continuous(limits = c(0,1), expand = c(0,0.01))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15))
  .plot10=.reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    geom_vline(xintercept = 0,size =1)+
    geom_hline(yintercept = 0,size =1)+
    geom_vline(xintercept = 0.05,colour="blue")+
    labs(y="Regression Coefficient",x="MAF (%)")+
    annotate("rect",xmin=0,xmax=1,ymin=.min,ymax=.max,alpha=0.2)+
    scale_x_continuous(limits = c(0,10), expand = c(0,0.15))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,0))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15))
  if(.in==1){
    .plot=cowplot::ggdraw()+
      cowplot::draw_plot(.plot10,x=0,y=0,width = 1,height = 1)+
      cowplot::draw_plot(.plot1+theme(axis.title = element_blank()),
                         x=0.45,y=0.15,width = 0.5,height = 0.5)
  }else if(.in==10){
    .plot=cowplot::ggdraw()+
      cowplot::draw_plot(.plot1,x=0,y=0,width = 1,height = 1)+
      cowplot::draw_plot(.plot10+theme(axis.title = element_blank()),
                         x=0.45,y=0.15,width = 0.5,height = 0.5)
  }
  plot(.plot)
  return(.plot)
}

}
