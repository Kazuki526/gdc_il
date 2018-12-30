# make rare variant focal patient list

# patient_list maked by prepare_tbl.R
patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")

all_patient_info=read_tsv("~/git/all_patient/all_patient_response.tsv",col_types = "cccccddccc") %>>%
  dplyr::rename(patient_id=submitter_id,age=diagnoses.0.age_at_diagnosis,gender=demographic.gender,
                race=demographic.race,ethnicity=demographic.ethnicity,stage=diagnoses.0.tumor_stage)

patient_race = all_patient_info %>>%
  mutate(race_=ifelse(race=="white" &ethnicity!="hispanic or latino","white",
                      ifelse(race=="black or african american" &ethnicity!="hispanic or latino",
                             "black","other"))) %>>%
  mutate(race_ = ifelse(is.na(race_),"other",race_)) %>>%
  dplyr::select(patient_id,race_) %>>%
  dplyr::rename(race=race_)

control_coverage = read_tsv("/Volumes/areca42TB/tcga/all_patient/by_patient_coverage_depth_control.tsv.gz") %>>%
  left_join(patient_list) %>>%
  filter(!is.na(gender))

#coverageどれくらいのcutoffが妥当？
ggplot(control_coverage)+geom_histogram(aes(x=called_bp),binwidth = 5000)+geom_vline(xintercept = c(600000,730000))

all_maf_for_cumulative_cont = read_tsv("/Volumes/areca42TB/tcga/all_patient/all_maf_for_cumulative_control.tsv.gz")
rare_variant_count = all_maf_for_cumulative_cont %>>%
  filter(ref!="-",alt!="-",MAF<0.0005)%>>%
  group_by(patient_id,mutype)%>>%
  summarise(variant_num=sum(MAC)) %>>%
  tidyr::spread(mutype,variant_num) %>>%
  mutate_all(funs(ifelse(is.na(.),0,.)))%>>%
  mutate(all_variant=missense+silent+splice+truncating)

  
gwas_patient = control_coverage %>>%
  filter(called_bp>600000,called_bp<730000) %>>%
  left_join(patient_race)%>>%
  dplyr::select(patient_id,cancer_type,called_bp,race)%>>%
  left_join(rare_variant_count) %>>%
  mutate_all(funs(ifelse(is.na(.),0,.)))%>>%
  mutate(variant_num_rank = cume_dist(all_variant))

write_df(gwas_patient,"/Volumes/areca42TB/tcga/all_patient/gwas_focal_patient_list.tsv")
