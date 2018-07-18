######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.R,quality_filter.Rを実行してからこのスクリプトを開始する！！！
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
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
####全cancer_typeまとめて####
cumulative_plot(.MAF_end = 0.5,.path = "all_cancer_type",.permu_file = "TSG/missense_05.tsv",
                .regression_size = 5,.pnum_size = 3)
#^^X.intercept =60.79655,R=-0.3837542,P=0
#cumulative_plot(.MAF_end = 0.23,.path = "all_cancer_type")
cumulative_plot(.MAF_end = 0.05,.path = "all_cancer_type",.permu_file = "TSG/missense_005.tsv",
                .regression_size = 5,.pnum_size = 3)
#cumulative_plot(.MAF_end = 1,.path = "all_cancer_type",.permu_file = "TSG/missense_1.tsv")
##^X.intercept =60.79254,R=-0.3124493,P=0
#cumulative_plot(.MAF_end = 5,.path = "all_cancer_type",.permu_file = "TSG/missense_5.tsv")
##^X.intercept =60.67493,R=-0.1491525,P=0
#cumulative_plot(.MAF_end = 10,.path = "all_cancer_type",.permu_file = "TSG/missense_10.tsv")
##^X.intercept =60.8105,R=-0.1417576,P=0
#cumulative_plot(.MAF_start = 0.23,.MAF_end = 0.5,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0,.MAF_end = 1,.path = "all_cancer_type")

#cumulative_plot(.MAF_start = 5,.MAF_end = 50,.path = "all_cancer_type", .more_5par = T)
#cumulative_plot(.MAF_start = 1,.MAF_end = 5,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.5,.MAF_end = 1,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.1,.MAF_end = 0.5,.path = "all_cancer_type")
#cumulative_plot(.MAF_start = 0.05,.MAF_end = 0.1,.path = "all_cancer_type")

#### cancer_type ごとに #####
.plot05_by  = cumulative_plot(.MAF_end = 0.5,.path = "by_cancer_type",.facet_by_cancer_type = T,
                              .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                              .permu_file = "TSG/missense_05_byCT.tsv")
#cumulative_plot(.MAF_end = 0.23,.path = "by_cancer_type",.facet_by_cancer_type = T)
.plot005_by =cumulative_plot(.MAF_end = 0.05,.path = "by_cancer_type",.facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                             .permu_file = "TSG/missense_005_byCT.tsv")
#cumulative_plot(.MAF_start = ,.MAF_end = ,.path = "by_cancer_type",.facet_by_cancer_type = T)

#########################################################################
# mutationのあるgene数でcount
cumulative_plot(.MAF_end = 0.05,.path = "gene",.by_gene = T,.permu_file = "TSG/check/missense_005_genenum.tsv")
#cumulative_plot(.MAF_end = 0.5,.path = "gene",.by_gene = T)
############################################################
##### 0.01%ごとのplot
#regression_table = make_regression_tabel(.max_maf = 50)
#write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn.tsv")
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside.pdf",.plot_reg,height = 6,width = 12)
# .plot_reg=regression_plot_in(regression_table,.in = 1)
# ggsave("age_plot/fig/regression/TSG_nonsyn_in1.pdf",.plot_reg,height = 6,width = 10)
# .plot_reg=regression_plot_in(regression_table,.in = 10)
# ggsave("age_plot/fig/regression/TSG_nonsyn_in10.pdf",.plot_reg,height = 6,width = 10)
if(0){
.plot_reg1=regression_table　%>>%
  regression_tbl_plot(.bl_ln = 0.05,.red_ln = 0.5)
ggsave("age_plot/cumulative/regression_plot-1.pdf",.plot_reg1,height = 6,width = 12)
.plot_reg10=regression_table　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.05,.red_ln = 0.5,.expand = 0.15)
ggsave("age_plot/cumulative/regression_plot-10.pdf",.plot_reg10,height = 6,width = 12)
}
##########################################################################################################
###################################### figrure用に調整 #######################################
.plot005=cumulative_plot(.MAF_end = 0.05,.path = "all_cancer_type",
                         .permu_file = "TSG/missense_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  ##MAF<1  : X.intercept =60.79254,R=-0.3124493,P=0
  geom_abline(aes(intercept=60.79254,slope=-0.3124493),colour="blue")+
  ##MAF<10 : X.intercept =60.8105,R=-0.1417576,P=0
  geom_abline(aes(intercept=60.8105,slope=-0.1417576),colour="green")
.plot005_by =cumulative_plot(.MAF_end = 0.05,.path = "by_cancer_type",.facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                             .permu_file = "TSG/missense_005_byCT.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.x = element_text(size =15),
                                      axis.title.y = element_text(size =15))+
                       ggtitle(label = NULL),
                     x=0,y=0.35,width=0.36,height=0.65)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_text(size = 15))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg.pdf",.plot,width = 14,height = 8)

# #MAF=0.5% (supply 用)
# .plot = cowplot::plot_grid(.plot05 + theme(axis.title = element_text(size =15),
#                                            title = element_text(size = 20)),
#                           .plot05_by + theme(axis.title.y = element_blank(),
#                                               axis.title.x = element_text(size = 15))+
#                             ggtitle(label = NULL),
#                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
#                           rel_widths = c(1,1.8))
# ggsave("age_plot/fig/TSG/maf05_nonsyn.pdf",.plot,width = 14,height = 8)
# .plot = cowplot::ggdraw()+
#   cowplot::draw_plot(.plot05 + theme(axis.title.x = element_text(size =15),
#                                      axis.title.y = element_text(size =15))+
#                        ggtitle(label = NULL),
#                      x=0,y=0.35,width=0.36,height=0.65)+
#   cowplot::draw_plot(.plot05_by + theme(axis.title.y = element_blank(),
#                                          axis.title.x = element_text(size = 15))+
#                        ggtitle(label = NULL),
#                      x=0.36,y=0,width=0.64,height=1)+
#   cowplot::draw_plot(regression_plot_byside(regression_table,.vline=0.5),
#                      x=0,y=0,width=0.36,height=0.35)+
#   cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
#   cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
#   cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
# 
# ggsave("age_plot/fig/TSG/maf05_nonsyn_withreg.pdf",.plot,width = 14,height = 8)
###regressionを下にしたfig2
# .plot = cowplot::ggdraw()+
#   cowplot::draw_plot(.plot005 + theme(axis.title.x = element_blank(),
#                                       axis.title.y = element_text(size =15))+
#                        ggtitle(label = NULL),
#                      x=0,y=0.4,width=0.36,height=0.6)+
#   cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
#                                          axis.title.x = element_text(size = 15))+
#                        ggtitle(label = NULL),
#                      x=0.36,y=0.4,width=0.64,height=0.6)+
#   cowplot::draw_text("Number of Nonsynonymous Variants",x=0.5,y=0.35,vjust = 0)+
#   cowplot::draw_plot(.plot_reg,x=0,y=0,width=1,height=0.35)+
#   cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
#   cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
#   cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
# 
# ggsave("age_plot/fig/TSG/maf005_nonsyn_withreg2.pdf",.plot,width = 14,height = 10)

#ここからは過去にボツったfigure
# .plot = cowplot::plot_grid(.plot05 + theme(axis.title = element_text(size =15),
#                                             title = element_text(size = 20)),
#                            .plot05_by + theme(axis.title.y = element_blank(),
#                                                axis.title.x = element_text(size = 15))+
#                              ggtitle(label = NULL),
#                            labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
#                            rel_widths = c(1,1.8))
# ggsave("age_plot/fig/TSG/maf05_nonsyn.pdf",.plot,width = 14,height = 8)
# .plot_reg = cowplot::ggdraw()+
#   cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+
#                        annotate("rect",xmin=0,xmax=1,ymin=-0.92,ymax=0.02,alpha=0.2)+
#                        scale_y_continuous(limits = c(-0.92,0.02),expand = c(0,0),
#                                           labels = c(-0.8,-0.6,-0.4,-0.2,0)),
#                      x=0.05, y=0.53, width = 0.9, height = 0.47)+
#   cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
#                      x=0.05, y=0   , width = 0.9, height = 0.53)+
#   cowplot::draw_text("regression coefficient",size = 15, x=0.025, y=0.5, angle=90)
# .plot_reg
# ggsave("age_plot/fig/presentation/TSG_nonsyn_reg.pdf",.plot_reg,width = 10,height =5 )
# 
# .plot05  = cumulative_plot(.MAF_end = 0.5,.save = F,.regression_size = 5,.pnum_size = 3,
#                            .permu_file = "TSG/missense_05.tsv")
# .plot005 = cumulative_plot(.MAF_end = 0.05,.save = F,.regression_size = 5,.pnum_size = 3,
#                            .permu_file = "TSG/missense_005.tsv")
# .plot = cowplot::plot_grid(.plot_reg,
#                    cowplot::plot_grid(.plot005 + theme(axis.title.x = element_blank(),
#                                                        axis.title.y = element_text(size = 15),
#                                                        axis.text = element_text(size = 12),
#                                                        title = element_text(size = 20)),
#                                       .plot05 + theme(title = element_text(size = 20),
#                                                       axis.title = element_text(size = 15),
#                                                       axis.text = element_text(size = 12)),
#                                       labels = c("b","c"),label_size = 30,
#                                       ncol = 1,rel_heights = c(1,1.1)),
#                    labels = c("a",""),ncol=2,rel_widths = c(1.5,1),label_size = 30)
# .plot
# ggsave("age_plot/fig/presentation/TSG_nonsyn_reg_and_violin.pdf",.plot,width = 15,height =7 )
################################################################################################
################################################################################################
################################################################################################
###################################### synonymou on TSG ########################################
#silent
# .plots05  = cumulative_plot(.MAF_end = 0.5,.mutype="silent",.path = "silent",
#                             .permu_file = "TSG/silent_05.tsv")
.plots005 = cumulative_plot(.MAF_end = 0.05,.mutype="silent",.path = "silent",
                            .permu_file = "TSG/silent_005.tsv")
.plots005_by = cumulative_plot(.MAF_end = 0.5,.path = "by_cancer_type",.facet_by_cancer_type = T,
                               .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.mutype = "silent",
                               .permu_file = "TSG/silent_005_byCT.tsv")
#cumulative_plot(.MAF_end = 0.01,.mutype="silent",.path = "silent")
############################################################
#regression_table_silent = make_regression_tabel(.mutype = "silent",.max_maf = 50)
#write_df(regression_table_silent,"age_plot/cumulative/regression/TSG_syn.tsv")
regression_table_silent = read_tsv("age_plot/cumulative/regression/TSG_syn.tsv")
.plot_regs = regression_table_silent %>>%
  regression_plot_byside(.red = NA)
##########################################################################################################
###################################### figrure用に調整 #######################################
.plot005=cumulative_plot(.MAF_end = 0.05,.path = "silent",.mutype = "silent",
                         .permu_file = "TSG/silent_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  ##MAF<1  : X.intercept =60.79254,R=-0.3124493,P=0
  geom_abline(aes(intercept=60.79254,slope=-0.3124493),colour="blue")+
  ##MAF<10 : X.intercept =60.8105,R=-0.1417576,P=0
  geom_abline(aes(intercept=60.8105,slope=-0.1417576),colour="green")
.plot005_by =cumulative_plot(.MAF_end = 0.05,.path = "silent",.mutype = "silent",.facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                             .permu_file = "TSG/silent_005_byCT.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.x = element_text(size =15),
                                      axis.title.y = element_text(size =15))+
                       ggtitle(label = NULL),
                     x=0,y=0.35,width=0.36,height=0.65)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_text(size = 15))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)

ggsave("age_plot/fig/TSG/maf005_syn_withreg.pdf",.plot,width = 14,height = 8)


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
  filter((role=="TSG"|role=="oncogene/TSG"),gene_symbol!="KMT2C") %>>%
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
cumulative_plot(.MAF_end = 0.06,.path = "onco",.role = "oncogene",
                .permu_file = "oncogene/missense05.tsv")
cumulative_plot(.MAF_end = 3,.path = "onco",.role = "oncogene",
                .permu_file = "oncogene/missense005.tsv")

cumulative_plot(.MAF_end = 0.02,.path = "onco/silent",.role = "oncogene",.mutype = "silent",
                .permu_file = )
####### missense by cancer type ########
cumulative_plot(.MAF_end = 0.06,.path = "onco/by_cancer_type",
                .role = "oncogene",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = )
cumulative_plot(.MAF_end = 3,.path = "onco/by_cancer_type",
                .role = "oncogene",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = )

##################################################################
regression_table_onco = make_regression_tabel(.role = "oncogene")
.plot_reg1=regression_table_onco　%>>%
  regression_tbl_plot(.bl_ln = 0.05)
ggsave("age_plot/cumulative/onco/regression_plot-1.pdf",.plot,height = 6,width = 12)
.plot_reg10=regression_table_onco　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.05,.red_ln = 2.5,.expand = 0.15)
ggsave("age_plot/cumulative/onco/regression_plot-10.pdf",.plot,height = 6,width = 12)
###################################### figrure用に調整 #######################################
.plot_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.098,ymax=0.28,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.098,0.28),expand = c(0,0)),
                     x=0.05, y=0.53, width = 0.9, height = 0.47)+
  cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
                     x=0.05, y=0   , width = 0.9, height = 0.53)+
  cowplot::draw_text("regression coefficient",size = 15, x=0.025, y=0.5, angle=90)

.plot25  = cumulative_plot(.MAF_end = 2.5, .role = "oncogene",
                            .save = F,.regression_size = 5,.pnum_size = 3)
.plot005 = cumulative_plot(.MAF_end = 0.05, .role = "oncogene",
                            .save = F,.regression_size = 5,.pnum_size = 3)
.plot = cowplot::plot_grid(.plot_reg,
                            cowplot::plot_grid(.plot005 + theme(axis.title.x = element_blank(),
                                                                 axis.title.y = element_text(size = 15),
                                                                 axis.text = element_text(size = 12),
                                                                 title = element_text(size = 20)),
                                               .plot25 + theme(title = element_text(size = 20),
                                                                axis.title = element_text(size = 15),
                                                                axis.text = element_text(size = 12)),
                                               labels = c("b","c"),label_size = 30,
                                               ncol = 1,rel_heights = c(1,1.1)),
                            labels = c("a",""),ncol=2,rel_widths = c(1.5,1),label_size = 30)
.plot
ggsave("age_plot/fig/presentation/oncogene_nonsyn_reg_and_violin.pdf",.plot,width = 15,height =7 )


###############################################################################
#silent
regression_table_silent_onco = make_regression_tabel(.mutype = "silent",.role = "oncogene")
.plots_reg1= regression_table_silent_onco %>>%
  regression_tbl_plot(.maf_max = 1, .bl_ln = 0.05, .red_ln = 0.02)
ggsave("age_plot/cumulative/onco/silent/regression_plot-1.pdf",.plot,height = 6,width = 12)
.plots_reg10= regression_table_silent_onco %>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.05,.red_ln = 0.02,.expand = 0.15)
ggsave("age_plot/cumulative/onco/silent/regression_plot-10.pdf",.plot,height = 6,width = 12)
###################################### figrure用に調整 #######################################
.plots_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plots_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.37,ymax=0.02,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.37,0.02),expand = c(0,0)),
                     x=0.05, y=0.53, width = 0.9, height = 0.47)+
  cowplot::draw_plot(.plots_reg1  + theme(axis.title.y = element_blank())+
                       scale_y_continuous(breaks = c(-0.3,-0.2,-0.1),labels = c(-0.3,-0.2,-0.1)),
                     x=0.05, y=0   , width = 0.9, height = 0.53)+
  cowplot::draw_text("regression coefficient",size = 15, x=0.025, y=0.5, angle=90)

.plots002  = cumulative_plot(.MAF_end = 0.02, .mutype = "silent",.role = "oncogene",
                            .save = F,.regression_size = 5,.pnum_size = 3)
.plots005 = cumulative_plot(.MAF_end = 0.05, .mutype = "silent",.role = "oncogene",
                            .save = F,.regression_size = 5,.pnum_size = 3)
.plots = cowplot::plot_grid(.plots_reg,
                            cowplot::plot_grid(.plots002 + theme(axis.title.x = element_blank(),
                                                                 axis.title.y = element_text(size = 15),
                                                                 axis.text = element_text(size = 12),
                                                                 title = element_text(size = 20)),
                                               .plots005 + theme(title = element_text(size = 20),
                                                                axis.title = element_text(size = 15),
                                                                axis.text = element_text(size = 12)),
                                               labels = c("b","c"),label_size = 30,
                                               ncol = 1,rel_heights = c(1,1.1)),
                            labels = c("a",""),ncol=2,rel_widths = c(1.5,1),label_size = 30)
.plots
ggsave("age_plot/fig/presentation/oncogene_syn_reg_and_violin.pdf",.plots,width = 15,height =7 )


##################################################################################################################
library(perm)
#MAF>5%で有意に発症年齢を下げているという結果が出たのでサイトごとに見てみる
onco_major=all_maf_for_cumlative%>>%
  filter(MAF>=0.05)%>>%
  left_join(truncating_count_onco) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
  filter((role=="oncogene"|role=="oncogene/TSG"),truncating_count_n==0) %>>%
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
  filter(role=="oncogene"|role=="oncogene/TSG") %>>%
  left_join(all_maf_for_cumlative %>>%filter(MAF<0.0006)%>>%
              group_by(gene_symbol,patient_id) %>>% summarise(missense_num=sum(MAC))) %>>%
  group_by(gene_symbol)%>>%mutate(focal=ifelse(length(patient_id)>1,"yes","no")) %>>%ungroup()%>>%
  nest(-gene_symbol,-role) %>>%
  mutate(regression = purrr::map(data, ~get_lm_coef(.))) %>>%
  dplyr::select(-role,-data) %>>%unnest() %>>%
  arrange(missense_num) %>>%
  write_df("~/Dropbox/install/tvz/oncogene_by_gene_regression_coef.tsv")
