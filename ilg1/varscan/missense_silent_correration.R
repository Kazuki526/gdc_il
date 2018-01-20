###nonsynonymou(missense)とsilent(synonymous)の数に相関があるのか？など
MAF=0.0001

missense_silent_num = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
              dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
  filter((mutype=="missense"|mutype=="silent") & MAF<MAF) %>>%
  group_by(cancer_type,patient_id,mutype) %>>%
  summarise(mutation_count=n()) %>>%ungroup()%>>%
  tidyr::spread(mutype,mutation_count) %>>%
  {left_join(patient_list,.)}%>>%
  mutate(missense=ifelse(is.na(missense),0,missense),
         silent  =ifelse(is.na(silent)  ,0,silent)) 
#ただただ相関をみてみる
missense_silent_num%>>%
{cor.test(.$missense,.$silent)}

#violin plotで
missense_silent_num %>>%
  ggplot(aes(x=as.factor(missense),y=silent))+
  geom_violin()+
  geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
#cancer type ごとに
missense_silent_num %>>%
  ggplot(aes(x=as.factor(missense),y=silent))+
  geom_violin()+
  geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2) +facet_wrap(~cancer_type)

#こっちのがみやすい？(legendをlogに治そう)
missense_silent_num %>>%
  count(missense,silent) %>>%
  mutate(log_n=log10(n)) %>>%
  ggplot(aes(x=missense,y=silent))+
  geom_raster(aes(fill=log_n))+scale_fill_distiller(palette = "BuGn",direction = 1)
#cancer type ごと
missense_silet_mean = missense_silent_num %>>%
  group_by(cancer_type) %>>%
  summarise(missense=mean(missense),silent=mean(silent))
missense_silent_num %>>%
  #count(missense,silent) %>>%
  count(cancer_type,missense,silent) %>>%
  mutate(log_n=log10(n)) %>>%
  ggplot(aes(x=missense,y=silent))+
  geom_raster(aes(fill=log_n))+scale_fill_distiller(palette = "BuGn",direction = 1)+
  geom_point(data=missense_silet_mean)+
  facet_wrap(~ cancer_type)


#変異の数はポアソン分布？
all_patient_num=6418
lambda_silent = missense_silent_num %>>%
  count(silent) %>>%
  {sum(.$n * .$silent)/all_patient_num}
lambda_missense=missense_silent_num %>>%
  count(missense) %>>%
  {sum(.$n * .$missense)/all_patient_num}
#missense
missense_silent_num %>>%
  count(missense) %>>%
  mutate(pois_expect=all_patient_num * dpois(missense,lambda_missense)) %>>%
  rename(observed=n) %>>%
  tidyr::gather(key,value,observed,pois_expect) %>>%
  ggplot(aes(x=missense,y=value))+
  geom_bar(aes(fill=key),stat = "identity",position = "dodge")
#silent
missense_silent_num %>>%
  count(silent) %>>%
  mutate(pois_expect=all_patient_num * dpois(silent,lambda_silent)) %>>%
  rename(observed=n) %>>%
  tidyr::gather(key,value,observed,pois_expect) %>>%
  ggplot(aes(x=silent,y=value))+
  geom_bar(aes(fill=key),stat = "identity",position = "dodge")

#ポアソン分布と仮定すると？
tibble(missense=rpois(all_patient_num,lambda_missense),silent=rpois(all_patient_num,lambda_silent)) %>>%
  count(missense,silent) %>>%
  {cor.test(.$missense,.$silent)}
tibble(missense=rpois(all_patient_num,lambda_missense),silent=rpois(all_patient_num,lambda_silent)) %>>%
  count(missense,silent) %>>%
  mutate(log_n=log10(n)) %>>%
  #count(cancer_type,missense,silent) %>>%
  ggplot(aes(x=missense,y=silent))+
  geom_raster(aes(fill=log_n))+scale_fill_distiller(palette = "BuGn",direction = 1)


tibble(all_count=rpois(all_patient_num,lambda_missense+lambda_silent)) %>>%
  mutate(missense=rbinom(all_patient_num,all_count,lambda_missense/(lambda_missense+lambda_silent))) %>>%
  mutate(silent=all_count - missense) %>>%
  count(missense,silent) %>>%
  {cor.test(.$missense,.$silent)}
tibble(all_count=rpois(all_patient_num,lambda_missense+lambda_silent)) %>>%
  mutate(missense=rbinom(all_patient_num,all_count,lambda_missense/(lambda_missense+lambda_silent))) %>>%
  mutate(silent=all_count - missense) %>>%
  count(missense,silent) %>>%
  mutate(log_n=log10(n)) %>>%
  ggplot(aes(x=missense,y=silent))+
  geom_raster(aes(fill=log_n))+scale_fill_distiller(palette = "BuGn",direction = 1)
