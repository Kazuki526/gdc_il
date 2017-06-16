library(tidyr)
library(plyr)
library(dplyr)
library(pipeR)
library(stringr)
library(ggplot2)
library(gridExtra)
library(readr)
library(readxl)
library(XML)
library(gtools)
library(purrr)

setwd("/Volumes/areca42TB/tcga/CNA/all_patient/")

############とりあえずもう使わない？
if(0){
cna_leng_data=read_tsv("all_patient_cna_length.tsv") %>>%
  mutate(age=round((diagnoses_age/365.25)*10)/10)

#とりあえずwhole genome duplicateは発症年齢と相関ありそう？
cna_leng_data %>>%
  ggplot()+
  geom_boxplot(aes(x=as.factor(WGD),y=age))

cna_leng_data %>>%
  #filter(line_num<100) %>>%
  filter(WGD==2) %>>%
  mutate(leng_del=leng_0 + leng_1 + leng_2,leng_amp=leng_2 + leng_3 + leng_4 + leng_5_ ) %>>%
  tidyr::gather(key = chr_n, value = leng, leng_del,leng_amp) %>>%
  ggplot(aes(x=leng,y=age))+
  geom_point()+
  facet_grid(chr_n ~ cancer_type)+
  stat_smooth(method = "lm", se = FALSE, colour = "red", size = 1)
}
#####################
cna_leng_data=read_tsv("all_patient_cna_length2.tsv") %>>%
  mutate(age=round((diagnoses_age/365.25)*10)/10)


