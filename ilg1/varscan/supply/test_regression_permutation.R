all_maf_for_cumulative = read_tsv("all_patient/all_maf_for_cumulative.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
patient_with_ps = read_tsv("/Volumes/areca42TB/tcga/all_patient/pathogenic_site_list.tsv")
snp_permu = function(.tbl,.times,.print=T){
  if(.times %% 1000 == 0 & .print) {print(paste0("permutation ",.times," times now"))}
  .tbl_sample=.tbl%>>%mutate(age=sample(age,length(age)))
  as.tibble(as.list(coef(lm(age/365.25 ~ missense_num,data = .tbl_sample))))
}


###############################################################################
permute_make_regression_tabel = function(.maf=all_maf_for_cumulative,.role = "TSG",
                                         .mutype="missense",.max_maf=10){
  regression_out = function(.class,.maf,.patient_list){
    ##missense の数
    if((.class*10000) %% 100 == 0){print(paste0("doing MAF=",.class*100))}
    missense_count = .maf %>>%
      filter(MAF <= .class)%>>%
      group_by(cancer_type,patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>% ungroup() %>>%
      {left_join(.patient_list,.,by = c("cancer_type","patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num)) 
    #相関直線を
    tibble::tibble(times=seq(1,10000,by=1)) %>>%
      mutate(tbl=purrr::map(times,~snp_permu(.tbl = missense_count,.,.print = F))) %>>%
      unnest() %>>%
      rename(regression = missense_num)
  }
  .patient_list=.maf %>>%
    left_join(driver_genes %>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter((mutype=="truncating"|mutype=="splice"),role==.role) %>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
    left_join(patient_with_ps,by=c("cancer_type","patient_id")) %>>%filter(is.na(significance))%>>%
    dplyr::select(-significance,-truncating_count_n)
  .maf = .maf %>>%filter(mutype==.mutype)%>>%
    left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
                dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
    filter(role==.role)
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.patient_list)))%>>%
    unnest()
}
permute_make_regression_tabel_cont = 
  function(.maf=all_maf_for_cumulative_cont,.role = "TSG",.mutype="missense",.max_maf=10){
    regression_out = function(.class,.maf){
      ##missense の数
      if((.class*10000) %% 100 == 0){print(paste0("doing MAF=",.class*100))}
      missense_count = .maf %>>%
        filter(MAF <= .class)%>>%
        group_by(cancer_type,patient_id) %>>%
        summarise(missense_num=sum(MAC)) %>>% ungroup() %>>%
        {left_join(patient_list,.,by = c("cancer_type","patient_id"))} %>>%
        mutate(missense_num = ifelse(is.na(missense_num),0,missense_num)) 
      #相関直線を
      tibble::tibble(times=seq(1,10000,by=1)) %>>%
        mutate(tbl=purrr::map(times,~snp_permu(.tbl = missense_count,.,.print = F))) %>>%
        unnest() %>>%
        rename(regression = missense_num)
    }
    .truncating_gene = all_maf_for_cumulative_cont %>>%
      filter(mutype=="truncating"|mutype=="splice") %>>%
      group_by(patient_id,gene_symbol) %>>%
      summarise(truncating_focal = "truncate") %>>%
      ungroup()
    .maf = .maf %>>%filter(mutype==.mutype)%>>%
      left_join(.truncating_gene,by = c("patient_id","gene_symbol")) %>>%
      filter(is.na(truncating_focal)) %>>%dplyr::select(-truncating_focal)
    tibble::tibble(MAF=1:(.max_maf*100)) %>>%
      mutate(MAF = MAF/10000) %>>%
      mutate(regression = purrr::map(MAF,~regression_out(.,.maf)))%>>%
      unnest()
}

perm_regression_nonsyn = permute_make_regression_tabel()
write_df(perm_regression_nonsyn,"age_plot/permute_tbl/by_maf/CDG_nonsyn.tsv.gz")

perm_regression_syn = permute_make_regression_tabel(.mutype = "silent")
write_df(perm_regression_syn,"age_plot/permute_tbl/by_maf/CDG_syn.tsv.gz")

perm_regression_nonsyn_cont = permute_make_regression_tabel()
write_df(perm_regression_nonsyn_cont,"age_plot/permute_tbl/by_maf/control_nonsyn.tsv.gz")

perm_regression_syn_cont = permute_make_regression_tabel(.mutype = "silent")
write_df(perm_regression_nonsyn_cont,"age_plot/permute_tbl/by_maf/control_syn.tsv.gz")


