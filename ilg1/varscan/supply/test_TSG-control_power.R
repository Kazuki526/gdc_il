loadNamespace('cowplot')
library(viridis)
all_maf_for_cumulative = read_tsv("all_patient/all_maf_for_cumulative.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("all_patient/all_maf_for_cumulative_control.tsv.gz")
patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
all_patient_info = read_tsv("~/git/all_patient/all_patient_response.tsv")
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  mutate(role=ifelse(role=="oncogene/TSG","TSG",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))
###############################################################################################################
.MAF=0.0005
#cancer gene(TSG)におけるnonsynonymouとsynonymousに相関？
missense_silent_num = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene,role)%>>%
              dplyr::rename(gene_symbol=gene), by = "gene_symbol") %>>%
  filter(role == "TSG") %>>%
  filter((mutype=="missense"|mutype=="silent") & MAF<.MAF) %>>%
  group_by(cancer_type,patient_id,mutype) %>>%
  summarise(mutation_count=sum(MAC)) %>>%ungroup()%>>%
  tidyr::spread(mutype,mutation_count) %>>%
  {left_join(patient_list,.)}%>>%
  mutate(missense=ifelse(is.na(missense),0,missense),
         silent  =ifelse(is.na(silent)  ,0,silent))
#controlのmissenseとcancer geneのmissense
missense_silent_num_cont = all_maf_for_cumulative_cont %>>%
  filter((mutype=="missense"|mutype=="silent") & MAF<.MAF) %>>%
  group_by(cancer_type,patient_id,mutype) %>>%
  summarise(mutation_count=sum(MAC)) %>>%ungroup()%>>%
  tidyr::spread(mutype,mutation_count) %>>%
  {left_join(patient_list,.)}%>>%
  mutate(missense_cont=ifelse(is.na(missense),0,missense),
         silent_cont  =ifelse(is.na(silent)  ,0,silent))%>>%
  dplyr::select(-missense,-silent)
##############################################################################################
#TSGにおけるmissense(nonsynonymous)とsilent(synonymous)の相関を見てみる
.plot = missense_silent_num %>>%
  ggplot(aes(x=missense,y=silent))+
   geom_violin(aes(group=missense),scale = "count",bw=0.5)+
   #geom_boxplot(width=.3,fill="black")+ 
   stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  stat_smooth(method=lm, formula = y ~ 0 +x,se=F,colour="black")+
  #stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
  xlab("Number of Nonsynonumous Variants")+ylab("Number of Synonymous Variants")+
  #ggtitle("TSG synonymous and nonsynonymous correlation")+
  theme_bw()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12))
.plot
ggsave("age_plot/fig/further_research/TSG_ns_corre.pdf",.plot,width = 5,height = 5)
#相関検定
cor.test(missense_silent_num$missense,missense_silent_num$silent,method = "pearson")
#P_value < 4.752e-13
#相関強

#ageをbinにしてnonsynonymous/synonymousの分布どうなっている？
.plot_all = missense_silent_num %>>%
  mutate(age = signif(age/365.25,4)) %>>%
  mutate(age_group = (age%/%5)*5 ) %>>%
  group_by(age_group) %>>%
  summarise(missense_sum =sum(missense),
            silent_sum=sum(silent),
            patient_num = n()) %>>%ungroup()%>>%
  filter(missense_sum>0,silent_sum>0) %>>%
  mutate(`nonsynonymous/synonymous` = missense_sum/silent_sum) %>>%
  ggplot(aes(x=reorder(as.factor(age_group),age_group),y=`nonsynonymous/synonymous`))+
  geom_point()+
  geom_text(aes(x=reorder(as.factor(age_group),age_group),y=0.2,label=patient_num),
            size=2, angle =-45)+
  theme_bw()+
  ggtitle("all cancer type")+xlab("age")
.plot_by = missense_silent_num %>>%
  mutate(age = signif(age/365.25,4)) %>>%
  mutate(age_group = (age%/%5)*5 ) %>>%
  group_by(cancer_type,age_group) %>>%
  summarise(missense_sum =sum(missense),
            silent_sum=sum(silent),
            patient_num = n()) %>>%ungroup()%>>%
  filter(missense_sum>0,silent_sum>0) %>>%
  mutate(`nonsynonymous/synonymous` = missense_sum/silent_sum) %>>%
  ggplot(aes(x=reorder(as.factor(age_group),age_group),y=`nonsynonymous/synonymous`))+
  geom_point()+
  geom_text(aes(x=reorder(as.factor(age_group),age_group),y=0.5,label=patient_num),
            size=2, angle =-45)+
  facet_wrap(~ cancer_type)+
  theme_bw()+xlab("age")+
  theme(axis.text.x = element_text(angle =-45,size=8))
.plot = cowplot::plot_grid(.plot_all, .plot_by, ncol = 2,labels = "auto", rel_widths = c(1,2) )
ggsave("age_plot/fig/further_research/TSG_nonsyn-syn.pdf",.plot,width = 13,height = 8)
##############################################################################################
#TSGにおけるmissense(nonsynonymous)とcontrolにおけるmissense(nonsynonymous)の相関を見てみる
.plot = missense_silent_num %>>%
  left_join(missense_silent_num_cont) %>>%
  ggplot(aes(x=missense,y=missense_cont))+
  geom_violin(aes(group=missense),scale = "count",bw=0.5)+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  stat_smooth(method=lm, formula = y ~ 0 +x,se=F,colour="black")+
  #stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
  xlab("Number of  TSG Nonsynonumous Variants")+ylab("Number of Normal gene Nonsynonymous Variants")+
  #ggtitle("TSG nonsynonymous and control gene nonsynonymous correlation")
  theme_bw()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12))
.plot
ggsave("age_plot/fig/further_research/nonsyn_TSG_control_corre.pdf",.plot,width = 5,height = 6)
#相関検定
tsg_cont_nonsyn = missense_silent_num %>>%
  left_join(missense_silent_num_cont)
cor.test(tsg_cont_nonsyn$missense,tsg_cont_nonsyn$missense_cont,method = "pearson")
#P_value < 2.2e-16
tsg_cont_nonsyn_lm = lm(missense ~ -1 + missense_cont, data = tsg_cont_nonsyn)
tsg_cont_nonsyn_lm$coefficients
1-pf(summary(tsg_cont_nonsyn_lm)$fstatistic["value"],summary(tsg_cont_nonsyn_lm)$fstatistic["numdf"],
     summary(tsg_cont_nonsyn_lm)$fstatistic["dendf"])
#相関強

#ageをbinにしてnonsynonymousのTSG/controlの分布どうなっている？
.plot_all = missense_silent_num %>>%
  left_join(missense_silent_num_cont) %>>%
  mutate(age = signif(age/365.25,4)) %>>%
  mutate(age_group = (age%/%5)*5 ) %>>%
  group_by(age_group) %>>%
  summarise(TSG_sum =sum(missense),
            control_sum=sum(missense_cont),
            patient_num = n()) %>>%ungroup()%>>%
  filter( TSG_sum>0,control_sum>0) %>>%
  mutate(`TSG/control_gene` = TSG_sum/control_sum) %>>%
  ggplot(aes(x=reorder(as.factor(age_group),age_group),y=`TSG/control_gene`))+
  geom_point()+
  geom_text(aes(x=reorder(as.factor(age_group),age_group),y=0.1,label=patient_num),
            size=2, angle = -45)+
  theme_bw()+
  xlab("age")+ggtitle("all cancer type")
.plot_by = missense_silent_num %>>%
  left_join(missense_silent_num_cont) %>>%
  mutate(age = signif(age/365.25,4)) %>>%
  mutate(age_group = (age%/%5)*5 ) %>>%
  group_by(cancer_type,age_group) %>>%
  summarise(TSG_sum =sum(missense),
            control_sum=sum(missense_cont),
            patient_num = n()) %>>%ungroup()%>>%
  filter( TSG_sum>0,control_sum>0) %>>%
  mutate(`TSG/control_gene` = TSG_sum/control_sum) %>>%
  ggplot(aes(x=reorder(as.factor(age_group),age_group),y=`TSG/control_gene`))+
  geom_point()+
  geom_text(aes(x=reorder(as.factor(age_group),age_group),y=0.1,label=patient_num),
            size=2, angle = -45)+
  facet_wrap(~ cancer_type)+
  theme_bw()+ xlab("age")+
  theme(axis.text.x = element_text(angle =-45,size=8))
.plot = cowplot::plot_grid(.plot_all, .plot_by, ncol = 2,labels = "auto", rel_widths = c(1,2) )
ggsave("age_plot/fig/further_research/nonsyn_TSG-control.pdf",.plot,width = 13,height = 10)
##################################################################################################
#tsg nonsynonymous とcontrol synonymous
cor.test(tsg_cont_nonsyn$missense,tsg_cont_nonsyn$silent_cont,method = "pearson")
tsg_nonsyn_cont_syn_lm =lm(silent_cont ~ -1 + missense ,data = tsg_cont_nonsyn)

tibble(to_=c("TSG Synonymous","Control Nonsynonymous","Control Synonymous"),
       R=c(tsg_nonsyn_syn_lm$coefficients,tsg_cont_nonsyn_lm$coefficients,tsg_nonsyn_cont_syn_lm$coefficients))
