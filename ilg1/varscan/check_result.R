#サプリにおいてresultがartifact出ないことを示すための図たち
################################################################
################# TSG のregression #############################
################################################################
##race = white
white_regression_table = make_regression_tabel(.race = "white")
.plot_regw = white_regression_table　%>>%
  regression_plot_byside()
.plot005w = cumulative_plot(.race = "white",.MAF_end = 0.05,.save = F,
                            .regression_size = 5,.pnum_size = 3,
                            .permu_file = "TSG/check/missense005_white.tsv")
.plot = cowplot::plot_grid(.plot_regw ,.plot005w + ggtitle(label = NULL),
                           labels = "auto",ncol = 1,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/white_regression.pdf",.plot,width = 8,height = 10)
##race = black
black_regression_table = make_regression_tabel(.race = "black")
.plot_regb = black_regression_table　%>>%
  regression_plot_byside()
.plot005b = cumulative_plot(.race = "black",.MAF_end = 0.05,.save = F,
                            .regression_size = 5,.pnum_size = 3,
                            .permu_file = "TSG/check/missense005_black.tsv")
.plot = cowplot::plot_grid(.plot_regw ,.plot005w + ggtitle(label = NULL),
                           labels = "auto",ncol = 1,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/white_regression.pdf",.plot,width = 8,height = 10)
.plot = cowplot::plot_grid(.plot_regb ,.plot005b + ggtitle(label = NULL),
                           labels = "auto",ncol = 1,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/black_regression.pdf",.plot,width = 8,height = 10)
#####これは使わない

##################################################################
### stage = 1,2でのみ
patient_stage12 = patient_list %>>%
  left_join(all_patient_info %>>%
              dplyr::select(patient_id,stage)) %>>%
  filter(stage==1)
regression_table_stage12 = make_regression_tabel(.patient_list = patient_stage12)
.plot_regw = regression_table_stage12　%>>%
  regression_plot_byside()
.plot005w = cumulative_plot(.MAF_end = 0.05,.save = F, .patient_list = patient_stage12,
                            .regression_size = 5,.pnum_size = 3,
                            .permu_file = "TSG/check/missense005_stage12.tsv")
.plot = cowplot::plot_grid(.plot_regw ,.plot005w + ggtitle(label = NULL),
                           labels = "auto",ncol = 1,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/stage12_regression.pdf",.plot,width = 8,height = 10)

##################################################################
### MAF = 0 をのぞいて見る
regression_table_remove0 = make_regression_tabel(.remove0 = T)
.plot_reg = regression_table_remove0　%>>%
  regression_plot_byside()
.plot005 = cumulative_plot(.MAF_start = 0.0001,.MAF_end = 0.05,.save = F,
                            .regression_size = 5,.pnum_size = 3,
                           .permu_file = "TSG/check/missense005_exclude0.tsv")
.plot = cowplot::plot_grid(.plot_reg ,.plot005 + ggtitle(label = NULL),
                           labels = "auto",ncol = 1,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/remove0_regression.pdf",.plot,width = 10,height = 10)

##################################################################
### vcfを1000genomeとUK10Kにしてみる
vcf_1kg = tally_1kg_all %>>%
  mutate(AF = ac_1kg/an_1kg) %>>%
  dplyr::select(chr,start,ref,alt,AF)
regression_table_1kg = make_regression_tabel(.vcf = vcf_1kg,.maf_filter = T)
.plot_reg_1kg = regression_table_1kg　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.05,.expand = 0.15)
.plot005_1kg = maf_trim_for_cumulative(.vcf = vcf_1kg)%>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3)
  
vcf_uk = uk10k %>>%
  mutate(AF = ac_uk/an_uk) %>>%
  dplyr::select(chr,start,ref,alt,AF)
regression_table_uk = make_regression_tabel(.vcf = vcf_uk,.maf_filter = T)
.plot_reg_uk = regression_table_uk　%>>%
  regression_tbl_plot(.maf_max = 10,.bl_ln = 0.05,.expand = 0.15)
.plot005_uk = maf_trim_for_cumulative(.vcf = vcf_uk) %>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3)

.plot = cowplot::plot_grid(.plot_reg_1kg + ggtitle("1000 genomes") +
                             theme(axis.title = element_text(size = 15),title = element_text(size = 20)),
                           .plot_reg_uk  + ggtitle("UK10K") +
                             theme(axis.title = element_text(size = 15),title = element_text(size = 20)),
                           .plot005_1kg + theme(title = element_text(size = 20),
                                                axis.title = element_text(size = 15),
                                                axis.text = element_text(size = 12)),
                           .plot005_uk + theme(title = element_text(size = 20),
                                                axis.title = element_text(size = 15),
                                                axis.text = element_text(size = 12)),
                           labels = "auto",ncol = 2,label_size = 30,rel_heights = c(1,1.5),
                           scale=0.9)
.plot
ggsave("age_plot/fig/check_result/regression_1kguk.pdf",.plot,width = 12,height = 8)

##################################################################
### quality filterしたものしていないもの
regression_table_dup = make_regression_tabel(.duplicate = F,.maf_filter = T)
.plot_reg_dup = regression_table_dup %>>% regression_plot_byside()
.plot005_dup = maf_trim_for_cumulative(.duplicate = F)%>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3,
                  .permu_file = "TSG/check/without_dup.tsv")

regression_table_som = make_regression_tabel(.somatic = F,.maf_filter = T)
.plot_reg_som = regression_table_som %>>% regression_plot_byside()
.plot005_som = maf_trim_for_cumulative(.somatic = F)%>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3,
                  .permu_file = "TSG/check/without_som.tsv")

regression_table_var = make_regression_tabel(.varscan = F,.maf_filter = T)
.plot_reg_var = regression_table_var　%>>%regression_plot_byside()
.plot005_var = maf_trim_for_cumulative(.varscan = F)%>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3,
                  .permu_file = "TSG/check/without_var.tsv")

regression_table_all = make_regression_tabel(.duplicate = F,.somatic = F,.varscan = F,.maf_filter = T)
.plot_reg_all = regression_table_all %>>% regression_plot_byside()
.plot005_all = maf_trim_for_cumulative(.duplicate = F,.somatic = F,.varscan = F)%>>%
  cumulative_plot(.MAF_end = 0.05,.save = F, .regression_size = 5,.pnum_size = 3,
                  .permu_file = "TSG/check/whithout_all.tsv")

.plot = cowplot::plot_grid(
  cowplot::plot_grid(.plot_reg_dup, .plot005_dup + ggtitle(label = NULL),
                     ncol = 2,rel_widths = c(2,1.3),scale=0.9),
  cowplot::plot_grid(.plot_reg_som, .plot005_som + ggtitle(label = NULL),
                     ncol = 2,rel_widths = c(2,1.3),scale=0.9),
  cowplot::plot_grid(.plot_reg_var, .plot005_var + ggtitle(label = NULL),
                     ncol = 2,rel_widths = c(2,1.3),scale=0.9),
  cowplot::plot_grid(.plot_reg_all, .plot005_all + ggtitle(label = NULL),
                     ncol = 2,rel_widths = c(2,1.3),scale=0.9),
  labels = "auto",label_size = 30,ncol = 1)
.plot
ggsave("age_plot/fig/check_result/qfilters_regression.pdf",.plot,width = 10,height = 16)


#ageをbinにして見て発症年齢低い人たちに引っ張られた回帰出ないことの証明
#MAF<0.05%
missense_count = all_maf_for_cumulative %>>%
  filter( MAF<=0.05/100,MAC!=0,mutype=="missense") %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role)%>>%dplyr::rename(gene_symbol=gene)) %>>%
  filter(role=="TSG") %>>%
  group_by(cancer_type,patient_id) %>>%
  summarise(missense_num=sum(MAC)) %>>%
  {left_join(patient_list,.)} %>>%
  mutate(missense_num = ifelse(is.na(missense_num),0,missense_num),
         age=round(age/365.25*100)/100)　%>>%
  left_join(truncating_count %>>%dplyr::select(-age)) %>>%
  filter(truncating_count_n==0) %>>%
  mutate(age_group_order = age %/% 5,
         age_group = paste0(as.character((age %/% 5)*5),"-",as.character((age %/% 5)*5+5)))
#all cancer type
missense_count %>>%
  ggplot(aes(x=reorder(age_group,age_group_order),y=missense_num))+
  geom_violin(scale = "count")+
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  geom_text(data = missense_count %>>% count(age_group),
            aes(x=age_group,y=-0.5,label = n),size = 3)+
  xlab("age onset")+ylab("number of misssense mutation")+
  theme_bw()
ggsave("age_plot/fig/check_result/missense_tsg_age_bin.pdf",width = 8,height = 5)
