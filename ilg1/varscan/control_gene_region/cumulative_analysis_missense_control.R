######！！！！！！！！注意！！！！！！！！########
# prepare_tbl_control_region.R を実行してからこのスクリプトを開始する！！！


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
             missense_count_order = missense_num, gender="male",
             age=-50, cancer_type="BRCA", patient_id="no_patient", truncating_count_n=0)
    return(rbind(.missense_count,.outcome))
  }else{return(.missense_count)}
}

######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
maf_trim_for_cumulative_cont = function(.maf=norm_maf_all_cont,.vcf=vcf_exac_cont,.race="all",.fdr=0.01,
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
    if(substitute(.vcf)=="vcf_exac_cont"){
      .vcf = .vcf %>>%
        mutate(AF=AC_Adj/AN_Adj)%>>%
        dplyr::select(chr,start,ref,alt,AF)
    }}else {stop(paste0(".race is wrong .race=",.race,"\nuse all, white or black!"))}
  ref_minor = ref_minor_focal_cont %>>%
    left_join(patient_race) %>>%
    filter(if(.race !="all"){race==.race}else{chr==chr}) %>>%
    left_join(left_join(tally_norm_maf_cont,.vcf)%>>%
                dplyr::select(-ac_cancer,-hom_cancer)) %>>%
    filter(AF > 0.5) %>>% mutate(MAF=1-AF) %>>%
    left_join(.maf) %>>%
    mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
           n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
           soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
           LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
                  mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,age,Protein_position) %>>%
    quality_filter_cont(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
    filter(MAC > 0)
  .maf %>>%
    quality_filter_cont(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    left_join(.vcf) %>>%
    mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
    filter(AF <= 0.5) %>>% mutate(MAF = AF) %>>%
    mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age,Protein_position) %>>%
    rbind(ref_minor) %>>%
    filter(!is.na(age),chr != "chrX")
}

all_maf_for_cumulative_cont = maf_trim_for_cumulative_cont()

#################################################################################################################
#0.01%ごとにregressionしてみる
make_regression_tabel = function(.maf=norm_maf_all_cont,.vcf=vcf_exac_cont,.race="all",
                                 .fdr=0.01,.mutype="missense",.max_maf=10,
                                 .database="all",.duplicate=T,.somatic=T,.varscan=T){
  regression_out = function(.minor_allele_frequency,.maf,.patient_list){
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .minor_allele_frequency)%>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>%
      {left_join(.patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num))
    #相関直線を
    lm=lm(age/365.25 ~ missense_num, data=missense_count)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"]))
  }
  if((.race=="all")&(substitute(.vcf)=="vcf_exac_cont")){
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
  if(substitute(.maf) =="norm_maf_all_cont"){
    .maf=maf_trim_for_cumulative_cont(.maf=.maf,.vcf=.vcf,.race=.race,.fdr=.fdr,.database=.database,
                                      .duplicate=.duplicate,.somatic=.somatic,.varscan=.varscan)
  }
  .maf = .maf %>>%filter(mutype==.mutype)#%>>%filter(MAF!=0)
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.patient_list)))%>>%
    unnest()
}
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

regression_table_cont = make_regression_tabel() 
.plot=regression_table_cont　%>>%
  regression_tbl_plot()
ggsave("age_plot/control_region/regression_plot-1_missense_control.pdf",.plot,height = 8,width = 20)
.plot=regression_table_cont　%>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_missense_control.pdf",.plot,height = 8,width = 20)

regression_table_silent_cont = make_regression_tabel(.mutype = "silent")
.plot= regression_table_silent_cont %>>%
  regression_tbl_plot(.maf_max = 1, .bl_ln = NULL, .red_ln = NULL)
ggsave("age_plot/control_region/regression_plot-1_silent_control.pdf",.plot,height = 8,width = 20)
.plot= regression_table_silent_cont %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_silent_control.pdf",.plot,height = 8,width = 20)

################################################
library(doParallel)
randome_sample_regression_missense = function(.repeat_num,.maf = all_maf_for_cumulative_cont,
                                              .mutype = "missense"){
  .tbl = .maf %>>% left_join(
      tibble(gene_symbol = sample(control_genes$gene_symbol,67)) %>>%
        mutate(focal = "ok")
    )%>>%
    filter(! is.na(focal),mutype == .mutype) %>>%
    make_regression_tabel(.max_maf = 1,.mutype = .mutype)
  if(.repeat_num %% 10 == 0){print(paste0("done ", .repeat_num, "times"))}
  return(.tbl)
}
repeat_regression_table_cont =
  tibble(repeat_num = 1:100) %>>%
  mutate(tbl = purrr::map(repeat_num, ~randome_sample_regression_missense(.))) %>>%
  unnest()
repeat_regression_table_cont_silent =
  tibble(repeat_num = 1:100) %>>%
  mutate(tbl = purrr::map(repeat_num, ~randome_sample_regression_missense(.mutype = "silent"))) %>>%
  unnest()
####################################################################################################
#violin plot 作ってみる
cumulative_plot = function(.maf=all_maf_for_cumulative_cont,.MAF_start = 0,.mutype="missense",
                           .MAF_end, .facet_by_cancer_type = F, .by_gene = F, .more_1par = F,
                           .race="all", .file_tail = ""){
  .patient_list=patient_list
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
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
    #missenseの数からmisssenseを持つ遺伝子の数に変えるか
    {if(.by_gene){.%>>%group_by(cancer_type,patient_id,gene_symbol) %>>%summarise(MAC=1)}else{.}} %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    #MAF>1%だと数が多いから3で割る
    {if(.more_1par){.%>>%mutate(missense_num = missense_num %/% 3)}else{.}} %>>%
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num)) %>>%
           {left_join(.patient_list,.)} %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order),
           missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
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
  write_df(regression,paste0("age_plot/control_region/violin_plot/lm_",.mutype,.MAF_start,"~",.MAF_end,"regression.tsv"))
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
  if(.more_1par){  #### MAF>1%の時
    .plot=.plot+scale_x_discrete(labels = legendx3)+
      xlab(paste0("number of ",.mutype," mutation"))
  }else{
    .plot = .plot+xlab(paste0("number of ",.mutype," mutation"))
  }
  .plot
  ggsave(paste0("age_plot/control_region/violin_plot/plot_",.mutype,.MAF_start,"~",.MAF_end,.file_tail,".pdf"),
         .plot,height = 16,width = 15)
}
cumulative_plot(.MAF_end = 0.05, .more_1par = T,.mutype = "silent")
cumulative_plot(.MAF_end = 0.05, .more_1par = T)
cumulative_plot(.MAF_end = 0.5, .more_1par = T,.mutype = "silent")
cumulative_plot(.MAF_end = 0.5, .more_1par = T)



################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
#可能性1:control geneにoncogeneが含まれていた。
#可能性2:essential geneに傷が入った状態でガンになる際passenger mutationがこれらの遺伝子で起こると
#        ガンですら生存が危うくなるため発症が遅くなる。
#可能性2を確かめるためhuman - mouse でKaKsが低い遺伝子のみで発症年齢をみてみよう！
enstXucsc = read_tsv("/Volumes/areca42TB/primates/ensembl2ucsc_Xref.tsv.gz",comment = "#")
primate_kaks = read_tsv("/Volumes/areca42TB/primates/ucid_primatesKAKS.tsv.gz")
control_genes_kaks = left_join(control_genes %>>%dplyr::select(gene_symbol,Transcript_ID),
                           enstXucsc,by = c("gene_symbol" = "HGNC_symbol")) %>>%
  filter(!is.na(UCSC_ID)) %>>%
  separate(UCSC_ID,into = c("ucsc_id","ucsc_id_num"),sep = "\\.") %>>%
  left_join(primate_kaks%>>%
              separate(ucid, into = c("ucsc_id","ucid_num"), sep = "\\.") %>>%
              dplyr::select(ucsc_id,ucid_num,Mouse)) %>>%
  filter(!is.na(Mouse)) %>>%
  separate(Mouse,into = c("kaks","cds_length"),sep = ":") %>>%
  mutate(kaks = parse_double(kaks), cds_length = parse_double(cds_length)) %>>%
  group_by(gene_symbol) %>>%
  summarise(cds_length = max(cds_length),kaks =kaks[which.max(cds_length)]) %>>%ungroup()

low_kaks_regression_table = all_maf_for_cumulative_cont %>>% 
  left_join(control_genes_kaks %>>%dplyr::select(gene_symbol,kaks,cds_length)) %>>%
  filter(!is.na(kaks),kaks <0.05,cds_length >100) %>>%
  make_regression_tabel()
.plot = low_kaks_regression_table %>>%
  regression_tbl_plot(.maf_max = 1,.bl_ln = NULL,.red_ln = NULL)
ggsave("age_plot/control_region/low_kaks/regression_plot-1_control.pdf",.plot,height = 8,width = 20)
.plot = low_kaks_regression_table %>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/low_kaks/regression_plot-10_control.pdf",.plot,height = 8,width = 20)

low_kaks_regression_table_silent = all_maf_for_cumulative_cont %>>% 
  left_join(control_genes_kaks %>>%dplyr::select(gene_symbol,kaks,cds_length)) %>>%
  filter(!is.na(kaks),kaks <0.05,cds_length >100 ) %>>%
  make_regression_tabel(.mutype = "silent")
.plot = low_kaks_regression_table_silent %>>%
  regression_tbl_plot(.maf_max = 1,.bl_ln = NULL,.red_ln = NULL)
ggsave("age_plot/control_region/low_kaks/regression_plot-1_silent_control.pdf",.plot,height = 8,width = 20)
.plot = low_kaks_regression_table_silent %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
ggsave("age_plot/control_region/low_kaks/regression_plot-10_silent_control.pdf",.plot,height = 8,width = 20)
