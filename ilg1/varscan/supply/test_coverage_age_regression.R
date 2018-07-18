control_genes = read_tsv("/Volumes/areca42TB/GRCh38_singlefasta/control_genes.tsv") %>>%
  filter(gene_symbol != "OR8U1") %>>% #GRCh37とGRCh38でゲノム上での向きが逆転しているから
  mutate(focal="yes")
control_genes_bed = read_tsv("/Volumes/areca42TB/GRCh38_singlefasta/control_genes_exon_with_splice_site.bed",
                             col_names = c("chr","start","end","gene_info","filter","strand"))
control_genes_exon_length = control_genes_bed %>>%
  mutate(length = end - start -5) %>>% {sum(.$length)}

cancergene_coverage = read_tsv("/Volumes/areca42TB/tcga/all_patient/by_patient_coverage_depth.tsv.gz") %>>%
  mutate(age_oncet = age/365.25) %>>%
  dplyr::rename(normal_mean_depth = norm_mean_depth) %>>%
  left_join(patient_list) %>>%
  filter(!is.na(gender))
control_coverage = read_tsv("/Volumes/areca42TB/tcga/all_patient/by_patient_coverage_depth_control.tsv.gz") %>>%
  mutate(normal_mean_depth = norm_mean_depth/control_genes_exon_length,
         tumor_mean_depth = tumor_mean_depth/control_genes_exon_length,
         age_oncet = age/365.25) %>>%
  left_join(patient_list) %>>%
  filter(!is.na(gender))
#####################################################################################################
#cancer driver gene
#TSGのdepthとageの相関
.lm_norm  = lm(age_oncet ~ normal_mean_depth, data = cancergene_coverage %>>%filter(role=="TSG"))
.lm_tumor = lm(age_oncet ~ tumor_mean_depth , data = cancergene_coverage %>>%filter(role=="TSG"))
regression_cancergene_tsg = 
  as.data.frame(as.list(coef(.lm_norm))) %>>%
  dplyr::rename(depth = normal_mean_depth) %>>%
  mutate(p_value_=(1 - pf(summary(.lm_norm)$fstatistic["value"],summary(.lm_norm)$fstatistic["numdf"],
                          summary(.lm_norm)$fstatistic["dendf"])),
         nt = "normal") %>>%
  bind_rows(as.data.frame(as.list(coef(.lm_tumor))) %>>%
              dplyr::rename(depth = tumor_mean_depth) %>>%
              mutate(p_value_=(1 - pf(summary(.lm_tumor)$fstatistic["value"],
                                      summary(.lm_tumor)$fstatistic["numdf"],
                                      summary(.lm_tumor)$fstatistic["dendf"])),
                     nt = "tumor"))
.plot = cancergene_coverage %>>%
  filter(role == "TSG") %>>%
  tidyr::gather(nt,depth,normal_mean_depth,tumor_mean_depth) %>>%
  mutate(nt=ifelse(nt == "normal_mean_depth","normal","tumor")) %>>%
  ggplot()+
  geom_point(aes(x=depth,y=age_oncet))+
  geom_abline(data = regression_cancergene_tsg,colour = "red",
              aes(intercept = X.Intercept.,slope = depth))+
  xlab("mean depth")+ylab("age onset")+
  facet_wrap( ~ nt) +
  theme_bw()
.plot
ggsave("age_plot/fig/check_result/TSG_depth_age_regression.pdf",.plot,width = 10,height = 5)

#oncogeneのdepthとageの相関
# .lm_norm = lm(age_oncet ~ normal_mean_depth, data = cancergene_coverage %>>%filter(role=="oncogene"))
# .lm_tumor= lm(age_oncet ~ tumor_mean_depth , data = cancergene_coverage %>>%filter(role=="oncogene"))
# regression_cancergene_onco = 
#   as.data.frame(as.list(coef(.lm_norm))) %>>%
#   dplyr::rename(depth = normal_mean_depth) %>>%
#   mutate(p_value_=(1 - pf(summary(.lm_norm)$fstatistic["value"],summary(.lm_norm)$fstatistic["numdf"],
#                           summary(.lm_norm)$fstatistic["dendf"])),
#          nt = "normal") %>>%
#   bind_rows(as.data.frame(as.list(coef(.lm_tumor))) %>>%
#               dplyr::rename(depth = tumor_mean_depth) %>>%
#               mutate(p_value_=(1 - pf(summary(.lm_tumor)$fstatistic["value"],
#                                       summary(.lm_tumor)$fstatistic["numdf"],
#                                       summary(.lm_tumor)$fstatistic["dendf"])),
#                      nt = "tumor"))
# .plot = cancergene_coverage %>>%
#   filter(role == "oncogene") %>>%
#   tidyr::gather(nt,depth,normal_mean_depth,tumor_mean_depth) %>>%
#   mutate(nt=ifelse(nt == "normal_mean_depth","normal","tumor")) %>>%
#   ggplot()+
#   geom_point(aes(x=depth,y=age_oncet))+
#   geom_abline(data = regression_cancergene_onco,colour = "red",
#               aes(intercept = X.Intercept.,slope = depth))+
#   xlab("mean depth")+ylab("age onset")+
#   facet_wrap( ~ nt) +
#   theme_bw()
# .plot
# ggsave("age_plot/fig/check_result/onco_depth_age_regression.pdf",.plot,width = 10,height = 5)
#########################################################################################################
#coverageとageの相関
.lm_tsg = lm(age_oncet ~ called_bp, data = cancergene_coverage %>>%filter(role=="TSG"))
regression_cancergene_coverage = 
  as.data.frame(as.list(coef(.lm_tsg))) %>>%
  mutate(p_value_=(1 - pf(summary(.lm_tsg)$fstatistic["value"],summary(.lm_tsg)$fstatistic["numdf"],
                          summary(.lm_tsg)$fstatistic["dendf"]))) 
.plot = cancergene_coverage %>>%
  filter(role == "TSG") %>>%
  ggplot()+
  geom_point(aes(x=called_bp,y=age_oncet))+
  geom_abline(data = regression_cancergene_coverage,colour = "black",
              aes(intercept = X.Intercept.,slope = called_bp))+
  xlab("coverage(bp)")+ylab("age onset")
.plot
ggsave("/Volumes/areca42TB/tcga/age_plot/fig/check_result/TSG_coverage_age_regression.pdf",.plot,
       width = 10,height = 5)

#mutation_num/coverageで相関を計算
tsg_missense_tbl = all_maf_for_cumulative %>>%
  filter(mutype == "missense") %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
              dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
  filter(role=="TSG" | role=="oncogene_TSG") %>>%
  left_join(truncating_count,by = c("cancer_type","patient_id")) %>>%
  filter(truncating_count_n==0)
regression_out = function(.minor_allele_frequency){
  ##missense の数
  missense_count = tsg_missense_tbl %>>%
    filter(MAF <= .minor_allele_frequency)%>>%
    group_by(cancer_type,patient_id) %>>%
    summarise(missense_num=sum(MAC)) %>>%
    {left_join(cancergene_coverage %>>% filter(role=="TSG"),.,
               by = c("patient_id", "cancer_type"))} %>>%
    mutate(missense_rate = ifelse(is.na(missense_num),0,missense_num/called_bp)) %>>%
    ungroup()
  #相関直線を
  lm=lm(age/365.25 ~ missense_rate, data=missense_count)
  as.data.frame(as.list(coef(lm))) %>>%
    mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                            summary(lm)$fstatistic["dendf"]))
}
regression_missense_rate = tibble::tibble(MAF=1:(10*100)) %>>%
  mutate(MAF = MAF/10000) %>>%
  mutate(regression = purrr::map(MAF,~regression_out(.)))%>>%
  unnest()
.plot
ggsave("age_plot/fig/check_result/regression_mutation_rate.pdf",width = 10,height = 5)
#####################################################################################################
#control gene region
.lm=lm(age_oncet ~ norm_mean_depth,data = control_coverage)
regression_control_ndepth=as.data.frame(as.list(coef(.lm))) %>>%
  mutate(p_value_=(1 - pf(summary(.lm)$fstatistic["value"],summary(.lm)$fstatistic["numdf"],
                          summary(.lm)$fstatistic["dendf"])))
cor.test(control_coverage$norm_mean_depth,control_coverage$age_oncet)
ggplot(data = control_coverage)+
  geom_point(aes(x=norm_mean_depth,y=age_oncet))+
  geom_abline(data = regression_control_ndepth,colour = "red",
              aes(intercept = X.Intercept.,slope = norm_mean_depth))

.lm=lm(age_oncet ~ tumor_mean_depth,data = control_coverage)
regression_control_tdepth=as.data.frame(as.list(coef(.lm))) %>>%
  mutate(p_value_=(1 - pf(summary(.lm)$fstatistic["value"],summary(.lm)$fstatistic["numdf"],
                          summary(.lm)$fstatistic["dendf"])))
ggplot(data = control_coverage)+
  geom_point(aes(x=tumor_mean_depth,y = age_oncet))

.lm=lm(age_oncet ~ called_bp,data = control_coverage)
regression_control_coverage=as.data.frame(as.list(coef(.lm))) %>>%
  mutate(p_value_=(1 - pf(summary(.lm)$fstatistic["value"],summary(.lm)$fstatistic["numdf"],
                          summary(.lm)$fstatistic["dendf"])))  

