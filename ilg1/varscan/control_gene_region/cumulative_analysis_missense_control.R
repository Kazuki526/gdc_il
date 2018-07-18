######！！！！！！！！注意！！！！！！！！########
# prepare_tbl_control_region.R を実行してからこのスクリプトを開始する！！！
loadNamespace('cowplot')

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
control_genes = read_tsv("/Volumes/areca42TB/GRCh38_singlefasta/control_genes.tsv") %>>%
  filter(gene_symbol != "OR8U1") %>>% #GRCh37とGRCh38でゲノム上での向きが逆転しているから
  mutate(focal="yes")

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
    filter(focal == "ok") %>>%dplyr::select(-focal)%>>%
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
    filter(MAC > 0)%>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,age,Protein_position)
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

#all_maf_for_cumulative_cont = maf_trim_for_cumulative_cont()
#write_df(all_maf_for_cumulative_cont,"all_patient/all_maf_for_cumulative_control.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
############################ ExAC data loading
vcf_exac_cont=read_tsv("/Volumes/areca42TB/exac/control_region/nontcga_liftovered_checked_likevcf.tsv.gz",
                       col_types = "cdccdddddddddddddddddddddddddddddddddd") %>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het))
vcf_exac_indel_cont=read_tsv("/Volumes/areca42TB/exac/control_region/nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                             col_types = "cdccdddddddddddddddddddddddddddddddddd")　%>>%
  mutate(AC_Het=ifelse(AC_Adj==AC_Het+AC_Hom*2,AC_Adj-AC_Hom*2,AC_Het)) 

vcf_exac_cont =rbind(vcf_exac_cont,vcf_exac_indel_cont) %>>%
  mutate(chr=paste0("chr",chr))
rm(vcf_exac_indel_cont)
#################################################################################################################
####################################################################################################
#violin plot 作ってみる
cumulative_plot_cont = function(.maf=all_maf_for_cumulative_cont,.MAF_start = 0,.mutype="missense",
                           .MAF_end, .facet_by_cancer_type = F, .by_gene = F, .more_1par = F,
                           .race="all", .file_tail = "",
                           .regression_size = 7,.pnum_size = 4,.save = T,
                           .permu=T,.test_tail="one",.permu_file=NA,.all_color="black"){
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
  .patient_list=patient_list
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
  ##missense の数
  missense_count = .maf %>>%
    filter(MAF>=.MAF_start/100, MAF<=.MAF_end/100,MAC!=0,mutype==.mutype) %>>%
    #missenseの数からmisssenseを持つ遺伝子の数に変えるか
    {if(.by_gene){.%>>%group_by(cancer_type,patient_id,gene_symbol) %>>%summarise(MAC=1)}else{.}} %>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    {left_join(.patient_list,.)} %>>%
    mutate(missense_num = ifelse(is.na(missense_num),0,missense_num),
           age=round(age/365.25*100)/100)
  #相関直線を
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
      regression = mutate(regression,p_value = .p)%>>%(?.)
    }else{
      regression %>>%
        mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                               summary(lm)$fstatistic["dendf"]))) %>>%(?.)%>>%
        mutate(cancer_type="All Cancer Types")
    }
  }
  if(.save){
    write_df(regression,paste0("age_plot/control_region/violin_plot/lm_",
                               .mutype,.MAF_start,"~",.MAF_end,"regression.tsv"))
  }
  legendx3=function(.legend){
    .legend = as.numeric(ifelse(.legend=="10-",10,.legend))
    .legend = .legend*3
    as.character(ifelse(.legend==30,"30-",.legend))
  }
  regression=regression %>>%
    mutate(out_reg=paste0("R=",signif(missense_num,2),", P",
                          ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  missense_count = missense_count %>>%
    #MAF>1%だと数が多いから3で割る
    {if(.more_1par){.%>>%mutate(missense_num = missense_num %/% 3)}else{.}} %>>%
    mutate(missense_count_n=as.character(ifelse(missense_num >= 10,"10-",missense_num)),
           missense_count_order=ifelse(missense_num >=10,10,missense_num)) %>>%
    mutate(missense_count_n = ifelse(is.na(missense_count_n),"0",missense_count_n),
           missense_count_order = ifelse(is.na(missense_count_order),0,missense_count_order))
  .max_count=max(missense_count$missense_count_order)
  .p_posi=.max_count/2 +1
  .coef_posi=.max_count+1
  ##バイオリンプロットで見やすく
  .plot = missense_count %>>%missense_filling()%>>%
    ggplot(aes(x=reorder(missense_count_n,missense_count_order), y=age))+
    geom_violin(scale = "count",colour=.all_color)+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill=.all_color,colour=.all_color,shape=21,size=2)+
    geom_abline(data = regression %>>%filter(p_value > 0.05),
                aes(intercept = X.Intercept.,slope = missense_num),linetype="dashed",colour=.all_color)+
    geom_abline(data = regression %>>%filter(p_value <= 0.05),
                aes(intercept = X.Intercept.,slope = missense_num),colour=.all_color)+
    ylim(0,90)+
    ggtitle(if(.MAF_start==0){paste0("MAF < ",.MAF_end," %")
    }else{paste0("MAF = ",.MAF_start," ~ ",.MAF_end," %")})+
    theme(title = element_text(size = 20), axis.title = element_text(size = 15),
          panel.grid.major.y = element_line(colour = "gray"),
          panel.grid.minor.y = element_line(colour = "gray"),
          panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.x = element_line(colour = "gray95"),panel.grid.minor.x = element_blank(),
          axis.ticks.y = element_blank())
  if(.facet_by_cancer_type){ #cancer_typeごとに
    .plot=.plot+
      facet_wrap( ~ cancer_type, ncol = 5)+
      geom_text(data = regression %>>%mutate(coefi=paste0("regression=",round(missense_num*100)/100)),
                aes(x=.coef_posi,y=15,label =coefi),size=3.5,hjust=1)+
      geom_text(data =missense_count %>>%count(cancer_type,missense_count_n),
                aes(x=missense_count_n,y=5,label=n),size=2.5,position="stack",angle=-45)+
      theme(axis.text = element_text(size=12), strip.text = element_text(size=12),
            strip.background = element_rect(fill="transparent", colour = "black"))
  }else{     #全cancer_type
    .plot=.plot+
      geom_text(data = regression,aes(x=.coef_posi,y=15,label =out_reg),
                 size=.regression_size,hjust=1,colour=.all_color)+
      geom_text(data =missense_count %>>%count(missense_count_n),colour=.all_color,
                aes(x=missense_count_n,y=5,label=n),size=.pnum_size,position="stack")+
      theme(axis.text = element_text(size=15))
  }
  .ns = ifelse(.mutype == "missense", "nonsynonymous","synonymous")
  if(.more_1par){  #### MAF>1%の時
    .plot=.plot+scale_x_discrete(labels = legendx3)+
      xlab(paste0("number of ",.ns," mutation"))
  }else{
    .plot = .plot+xlab(paste0("number of ",.ns," mutation"))
  }
  if(.save){
    ggsave(paste0("age_plot/control_region/violin_plot/plot_",
                  .mutype,.MAF_start,"~",.MAF_end,.file_tail,".pdf"),
           .plot,height = 8,width = 8)
  }
  plot(.plot)
  return(.plot)
}
###########################################################################################################
#0.01%ごとにregressionしてみる
make_regression_tabel_cont = function(.maf=all_maf_for_cumulative_cont,.vcf=vcf_exac_cont,.race="all",
                                 .fdr=0.01,.mutype="missense",.max_maf=50,.filter_maf=F,
                                 .database="all",.duplicate=T,.somatic=T,.varscan=T,
                                 .patient_list=patient_list){
  regression_out = function(.minor_allele_frequency,.maf,.patient_list){
    if((.minor_allele_frequency*10000) %% 100 == 0){print(paste0("doing MAF=",.minor_allele_frequency*100))}
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .minor_allele_frequency)%>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>%
      {left_join(.patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num)) %>>%
      ungroup()
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
  .patient_list=.patient_list
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      left_join(patient_race) %>>%
      filter(race==.race)
  }
  if(.filter_maf){
    .maf=maf_trim_for_cumulative_cont(.maf=.maf,.vcf=.vcf,.race=.race,.fdr=.fdr,.database=.database,
                                      .duplicate=.duplicate,.somatic=.somatic,.varscan=.varscan)
  }
  .truncating_gene = all_maf_for_cumulative_cont %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    group_by(patient_id,gene_symbol) %>>%
    summarise(truncating_focal = "truncate") %>>%
    ungroup()
  .maf = .maf %>>%filter(mutype==.mutype) %>>%
    left_join(.truncating_gene,by = c("patient_id","gene_symbol")) %>>%
    filter(is.na(truncating_focal)) %>>%dplyr::select(-truncating_focal)
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.patient_list)))%>>%
    unnest()
}

####################################################################################################
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",
                .permu_file = "control/silent005.tsv")
.plot005 = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,
                           .permu_file = "control/missense005.tsv")
cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,.mutype = "silent",
                .permu_file = "control/silent05.tsv")
.plot05 = cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,
                          .permu_file = "control/missense05.tsv")

#by cancer type
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.facet_by_cancer_type=T,
                .permu_file = "control/silent005_byCT.tsv")
.plot005_by = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.facet_by_cancer_type=T,
                              .permu_file = "control/missense005_byCT.tsv")
cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,.mutype = "silent",.facet_by_cancer_type=T,
                .permu_file = "control/silent05_byCT.tsv")
.plot05_by = cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,.facet_by_cancer_type=T,
                             .permu_file = "control/missense05_byCT.tsv")
####################################################################################################

# regression_table_cont = make_regression_tabel_cont() 
# write_df(regression_table_cont,"age_plot/control_region/regression_table/nonsyn.tsv")
# regression_table_silent_cont = make_regression_tabel_cont(.mutype = "silent")
# write_df(regression_table_silent_cont,"age_plot/control_region/regression_table/syn.tsv")

regression_table_cont=read_tsv("age_plot/control_region/regression_table/nonsyn.tsv")
regression_table_silent_cont=read_tsv("age_plot/control_region/regression_table/syn.tsv")
##############################################################################################
############################### figure 用に調整 ##############################################
.plot = cowplot::plot_grid(.plot005 + theme(axis.title = element_text(size =15),
                                            title = element_text(size = 20)),
                           .plot005_by + theme(axis.title.y = element_blank(),
                                               axis.title.x = element_text(size = 15))+
                             ggtitle(label = NULL),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
.plot
ggsave("age_plot/fig/control/maf005_nonsyn.pdf",.plot,width = 14,height = 8)

.plot = cowplot::plot_grid(.plot05 + theme(axis.title = element_text(size =15),
                                            title = element_text(size = 20)),
                           .plot05_by + theme(axis.title.y = element_blank(),
                                               axis.title.x = element_text(size = 15),
                                               axis.text.x = element_text(size=8))+
                             ggtitle(label = NULL),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
.plot
ggsave("age_plot/fig/control/maf05_nonsyn.pdf",.plot,width = 14,height = 8)




##################################################################
.plot_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.49,ymax=0.02,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.49,0.02),expand = c(0,0)),
                     x=0.05, y=0.53, width = 0.9, height = 0.45)+
  cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
                     x=0.05, y=0   , width = 0.9, height = 0.50)+
  cowplot::draw_text("regression coefficient",size = 30, x=0.025, y=0.5, angle=90)
.plot_reg
.plots_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plots_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.72,ymax=0.02,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.72,0.02),expand = c(0,0)),
                     x=0, y=0.53, width = 0.9, height = 0.45)+
  cowplot::draw_plot(.plots_reg1  + theme(axis.title.y = element_blank()),
                     x=0, y=0   , width = 0.9, height = 0.50)

.plot = cowplot::plot_grid(.plot_reg,.plots_reg)
ggsave("age_plot/fig/presentation/control_reg.pdf",.plot,width = 15,height = 8)

.plots005 = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",
                            .save = F,.regression_size = 5,.pnum_size = 3)+
  theme(title = element_text(size = 15),axis.text = element_text(size = 15),
        axis.title = element_text(size =20))
.plot005 = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,
                           .save = F,.regression_size = 5,.pnum_size = 3)+
  theme(title = element_text(size = 15),axis.text = element_text(size = 15),
        axis.title = element_text(size =20))
.plot = cowplot::plot_grid(.plot_reg,.plots_reg,.plot005,.plots005,
                   labels = "auto",label_size = 30,ncol = 2)
.plot
ggsave("age_plot/fig/control/ns_reg_and_violin.pdf",.plot,width = 15,height = 15)
###############################################
#stage1,2でのみ
regression_table_stage12_cont = make_regression_tabel_cont(.patient_list = patient_list %>>%
                                                        left_join(all_patient_info %>>%
                                                                    dplyr::select(patient_id,stage)) %>>%
                                                        filter(stage==1)) 
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot()
ggsave("age_plot/control_region/regression_plot-1_missense_stage1_control.pdf",.plot,height = 8,width = 20)
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_missense_stage_1_control.pdf",.plot,height = 8,width = 20)

regression_table_stage12_silent_cont = make_regression_tabel_cont(.patient_list = patient_list %>>%
                                                        left_join(all_patient_info %>>%
                                                                    dplyr::select(patient_id,stage)) %>>%
                                                        filter(stage==1), .mutype = "silent") 
.plot=regression_table_stage12_silent_cont　%>>%
  regression_tbl_plot()
ggsave("age_plot/control_region/regression_plot-1_silent_stage1_control.pdf",.plot,height = 8,width = 20)
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_silent_stage_1_control.pdf",.plot,height = 8,width = 20)


