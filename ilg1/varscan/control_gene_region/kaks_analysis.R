loadNamespace('cowplot')
control_genes = read_tsv("/Volumes/areca42TB/GRCh38_singlefasta/control_genes.tsv") %>>%
  filter(gene_symbol != "OR8U1") %>>% #GRCh37とGRCh38でゲノム上での向きが逆転しているから
  mutate(focal="yes")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
#可能性1:control geneにoncogeneが含まれていた。
#可能性2:essential geneに傷が入った状態でガンになる際passenger mutationがこれらの遺伝子で起こると
#        ガンですら生存が危うくなるため発症が遅くなる。
#可能性2を確かめるためhuman - mouse でKaKsが低い遺伝子のみで発症年齢をみてみよう！
mean_kaks = 0.193
# human - mouse の平均KaKsは0.192551524673489
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
  summarise(cds_length = max(cds_length),kaks =kaks[which.max(cds_length)]) %>>%ungroup() %>>%
  filter(cds_length>100)
.plot = control_genes_kaks %>>%
  ggplot()+
  geom_histogram(aes(x=kaks),binwidth = 0.01)+
  geom_vline(xintercept = mean_kaks,color="blue")+
  geom_hline(yintercept = 0)+
  theme_bw()
ggsave("age_plot/fig/further_research/control_gene_kaks.pdf",.plot,height = 5,width = 10)

# low_kaks_regression_table = all_maf_for_cumulative_cont %>>% 
#   left_join(control_genes_kaks %>>%dplyr::select(gene_symbol,kaks,cds_length)) %>>%
#   filter(!is.na(kaks),kaks < mean_kaks,cds_length >100) %>>%
#   make_regression_tabel_cont()
# high_kaks_regression_table = all_maf_for_cumulative_cont %>>% 
#   left_join(control_genes_kaks %>>%dplyr::select(gene_symbol,kaks,cds_length)) %>>%
#   filter(!is.na(kaks),kaks > mean_kaks,cds_length >100) %>>%
#   make_regression_tabelcont()
# 
# .plot_low = low_kaks_regression_table %>>%
#   regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
# .plot_high = high_kaks_regression_table %>>%
#   regression_tbl_plot(.maf_max = 10,.bl_ln = NULL,.red_ln = NULL,.expand = 0.15)
# 
# .plot = cowplot::plot_grid(.plot_low +ggtitle("low kaks control gene (271 gene)")+
#                              theme(axis.text = element_text(size = 10),axis.title = element_text(size = 15)),
#                            .plot_high +ggtitle("high kaks control gene (354 gene)")+
#                              theme(axis.text = element_text(size = 10),axis.title = element_text(size = 15)),
#                            ncol = 1, labels = "auto",label_size = 20)
# .plot
# ggsave("age_plot/fig/further_research/control_hl_kaks_regression.pdf",height = 10,width = 10)
############################################################################################################
maf_focal=all_maf_for_cumulative_cont%>>%
  filter(MAF<0.0005,mutype=="missense")
lm_reg_coef_cont=function(.gene_list){
  .missense_count = maf_focal%>>%inner_join(.gene_list,by="gene_symbol")%>>%
    group_by(patient_id) %>>%
    summarise(missense_num=sum(MAC))%>>%
    right_join(patient_list,by="patient_id")%>>%
    mutate(missense_num=ifelse(is.na(missense_num),0,missense_num))
  lm(age/365.25 ~ missense_num, data = .missense_count)$coefficients["missense_num"]
}
toplow_lm_r=function(.kaks_list,.gene_num=100){
  .top=.kaks_list %>>%
    mutate(rank=min_rank(kaks))%>>%
    filter(rank<=.gene_num)%>>%
    lm_reg_coef_cont()
  .low=.kaks_list %>>%
    mutate(rank=min_rank(desc(kaks)))%>>%
    filter(rank<=.gene_num)%>>%
    lm_reg_coef_cont()
  tibble::tibble(top=as.numeric(.top),low=as.numeric(.low))
}
#100
observe100=toplow_lm_r(control_genes_kaks)
permute_lm100 = tibble::tibble(times=1:10000) %>>%
  mutate(toplow = purrr::pmap(., function(times){
    if(times %%1000 ==0){print(paste0("now ",times," times"))}
    toplow_lm_r(control_genes_kaks%>>%
                mutate(kaks=sample(kaks,length(kaks))))
  }))%>>%unnest()
permute_lm100%>>%mutate(toplow=top/low)%>>%
  filter(toplow > observe100$top/observe100$low)
#200
observe200=toplow_lm_r(control_genes_kaks,200)
permute_lm200 = tibble::tibble(times=1:10000) %>>%
  mutate(toplow = purrr::pmap(., function(times){
    if(times %%1000 ==0){print(paste0("now ",times," times"))}
    toplow_lm_r(control_genes_kaks%>>%
                  mutate(kaks=sample(kaks,length(kaks))),200)
  }))%>>%unnest()
permute_lm200%>>%mutate(toplow=top-low)%>>%
  filter(toplow < observe200$top-observe200$low)
permute_lm200 %>>%ggplot()+
  geom_histogram(aes(x=top/low))

###############################################################################################################
plot_lm_sliding_kaks_rank = function(.w_size){
  gene_num=length(control_genes_kaks$gene_symbol)
  .R=tibble(group=c(0:(( gene_num-.w_size+10)%/%10))) %>>%
    mutate(min_rank=group*10 +1,
           max_rank=ifelse(group*10 +.w_size>gene_num,gene_num,group*10 +.w_size)) %>>%
    mutate(reg_coef = purrr::pmap_dbl(.,function(min_rank,max_rank,...){
      lm_reg_coef_cont(control_genes_kaks%>>%
                         mutate(rank=row_number(kaks))%>>%
                         filter(rank>=min_rank,rank<=max_rank))
    })) %>>%
    mutate(group_=paste0(min_rank,"-",max_rank)) %>>%
    ggplot()+
    geom_point(aes(x=reorder(group_,group),y=reg_coef))+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
}
.plot = plot_lm_sliding_kaks_rank(200)+
  labs(x= "Rank of Normal genes KaKs",y= "Regression Coeffisient")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("age_plot/fig/further_research/control_kaks_slidingwindow.pdf",width = 7,height = 4)    

plot_lm_sliding_kaks = function(.w_size){
  tibble(min_kaks=c(0:(1/.w_size))) %>>%
    mutate(max_kaks = min_kaks +.w_size) %>>%
    mutate(reg_coef=)
}
meanvar_by_patient=function(.num_tbl){
  patient_list%>>%left_join(.num_tbl,by="patient_id")%>>%
    mutate(n=ifelse(is.na(n),0,n))%>>%
    {tibble(mean=mean(.$n),variance=var(.$n))}
}
maf_focal %>>%count(gene_symbol,patient_id)%>>%
  nest(-gene_symbol)%>>%
  mutate(data_ = purrr::map(data,~meanvar_by_patient(.))) %>>%
  dplyr::select(-data)%>>%unnest() %>>%
  mutate(value=mean/variance)%>>%
  ggplot()+
  geom_point(aes(x=mean,y=variance))
  geom_histogram(aes(x=value))
  left_join(control_genes_kaks)%>>%mutate(rank=min_rank(kaks))%>>%View()
  gather(meanvar,value,-gene_symbol) %>>%
  ggplot()+
  geom_histogram(aes(x=value))+
  facet_grid( meanvar ~ . )

  .mutation=tibble(group=c(0:(( gene_num-.w_size+10)%/%10))) %>>%
    mutate(min_rank=group*10 +1,
           max_rank=ifelse(group*10 +.w_size>gene_num,gene_num,group*10 +.w_size)) %>>%
    mutate(mut_num = purrr::pmap_dbl(.,function(min_rank,max_rank,...){
      maf_focal %>>% inner_join(control_genes_kaks%>>%mutate(rank=min_rank(kaks))%>>%
                                  filter(rank >= min_rank, rank <= max_rank))
    })) %>>%
    mutate(group_=paste0(min_rank,"-",max_rank)) %>>%
    ggplot()+
    geom_line(aes(x=reorder(group_,group),y=mut_num))+
    theme(axis.text.x = element_text(angle = 90))
  cowplot::plot_grid(.R,.mutation,.ncol=1,rel_heights = c(1,1.2))