#this script is /working/ploting_1kg.R
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
setwd('/working/1000genomes/maf/')
write_df_watal = function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}
topdriver_bed=read_tsv("/Volumes/areca42TB/tcga/maf_norm/top_driver105.bed",col_names = c("chr","gene_start","gene_end","ids","score","strand"))%>>%
  tidyr::separate(ids,c("gene_symbol","ensg"),sep=";")
#brca_maf=read_tsv("/working/maf/7a07a833-4eab-44c9-bbf6-a64bd51c012e/TCGA.BRCA.muse.7a07a833-4eab-44c9-bbf6-a64bd51c012e.protected.maf.gz",comment = "#")
classify_consequence = function(.data) {
  dplyr::mutate(.data, mutype= dplyr::recode(Consequence,
                                             downstream_gene_variant = 'flank',
                                             `3_prime_UTR_variant` = 'flank',
                                             upstream_gene_variant = 'flank',
                                             `5_prime_UTR_variant` = 'flank',
                                             frameshift_variant = 'truncating',
                                             frameshift_variant = 'truncating',
                                             inframe_deletion = 'inframe_indel',
                                             inframe_insertion = 'inframe_indel',
                                             intron_variant = 'silent',
                                             splice_region_variant = 'splice',
                                             coding_sequence_variant = 'missense',
                                             missense_variant = 'missense',
                                             stop_gained = 'truncating',
                                             stop_lost = 'truncating',
                                             stop_retained_variant = 'silent',
                                             synonymous_variant = 'silent',
                                             splice_acceptor_variant = 'splice',
                                             splice_donor_variant = 'splice',
                                             protein_altering_variant = 'missense',
                                             start_lost = 'truncating',
                                             `splice_region_variant,intron_variant` = 'splice',
                                             `stop_gained,frameshift_variant` = 'truncating',
                                             `splice_region_variant,synonymous_variant`='splice',
                                             `splice_region_variant,5_prime_UTR_variant`='splice',
                                             `missense_variant,splice_region_variant`='missense',
                                             `intron_variant,non_coding_transcript_variant`='silent',
                                             `non_coding_transcript_exon_variant,non_coding_transcript_variant`='silent'))
}
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%
  mutate(role=ifelse(role=="oncogene/TSG","oncogene",role)) %>>%
  mutate(role=ifelse(is.na(role),"TSG",role))
plot_posi=read_tsv("~/ascat/data/plot_posi.tsv")

maf_root='/working/1000genomes/maf'
.infiles = list.files(maf_root, '.maf$', recursive=TRUE, full.names=TRUE) %>>% (?.)

maf_list = data_frame(path=.infiles, base=basename(path)) %>>%
  tidyr::separate(base, c('sample_allele', 'no'), '\\.', extra='drop') %>>%
  dplyr::select(-no) %>>%
  tidyr::separate(sample_allele,c("sample","allele"),'_',extra = 'drop') %>>%
  (?.)

.colnames = c('Hugo_Symbol','Chromosome','Start_Position','Consequence','PolyPhen')
.cols = .colnames %>>%
{setNames(c('c','c','d','c','c'), .)} %>>%
{do.call(readr::cols_only, as.list(.))}

strip_maf = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    filter(!mutype=='silent') %>>%
    filter(!mutype=='flank')
}

if(0){#test
  .colnamestest = c('Hugo_Symbol','Chromosome','Start_Position','Consequence','Protein_position')
  .colstest = .colnamestest %>>%
  {setNames(rep('c', length(.)), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  all_strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.colstest) %>>%
      classify_consequence() %>>%
      tidyr::separate(Protein_position,c("pposi","pall"),"/")
  }
  mafstest=maf_list %>>%
    mutate(maf=purrr::map(path,~all_strip_maf(.))) %>>%
    select(-path) %>>%
    unnest()
  mafstest %>>%
    filter(!pposi=="-") %>>%
    group_by(sample,allele,Hugo_Symbol,pposi) %>>%
    tally ()%>>%
    filter(n>1) %>>%View() #how many?
  
} 


mafs=maf_list %>>% #head(105) %>>%  
  mutate(maf=purrr::map(path,~strip_maf(.))) %>>%
  select(-path) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  mutate(gene_width=gene_end - gene_start) %>>%
  filter(start>gene_start,start<gene_end) %>>%
  mutate(freqposi=(start - gene_start)/gene_width) %>>%
  mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi)) %>>%
  dplyr::select(sample,allele,gene_symbol,freqposi,mutype) %>>%
  left_join(plot_posi) %>>%
  mutate(freqposi=freqposi + gposi) %>>%
  dplyr::select(-gposi) %>>%
  mutate(allele=ifelse(allele=="allele1",0,1))

#.maf=mafs%>>%filter(sample=='HG00096',groupe=="chr1:3")%>>%group_by(sample,groupe)%>>%nest(.key=variations)
plot_by_sample=function(.maf){
  sample=.maf$sample
  chrgroupe=.maf$groupe
  .plt=.maf$variations[[1]] %>>%
    mutate(quoti=(freqposi %/% 0.05)+0.025,surplus= freqposi %% 0.05) %>>%
    group_by(allele,quoti) %>>%
    mutate(rank=dplyr::row_number(surplus)) %>>%
    ungroup() %>>%
    mutate(allele = allele+((rank-1)*0.12+0.04),freqposi=quoti / 20 +0.025) %>>%
    dplyr::select(-quoti,-surplus) %>>%
    ggplot()+
    geom_rect(data = plot_posi%>>%filter(groupe==chrgroupe),aes(xmin=gposi,xmax=gposi+1,ymin=0,ymax=1),fill="gray")+
    geom_rect(data = plot_posi%>>%filter(groupe==chrgroupe),aes(xmin=gposi,xmax=gposi+1,ymin=1,ymax=2),fill="gray50")+
    facet_grid(groupe~.)+
    coord_cartesian(ylim = c(-0.4,2.2),xlim=c(0,21),expand=F)+
    geom_hline(yintercept = c(1,2),colour="white")+
    geom_vline(xintercept = 0:20,colour="black",size=0.1)+
    annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gene_symbol},
             x=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gposi} + 0.5,
             y=-0.15, size = 1.5,colour="red")+
    annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gene_symbol},
             x=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gposi} + 0.5,
             y=-0.15, size = 1.5,colour="black")+
    geom_point(aes(x=freqposi,y= allele,fill=mutype),size=1,stroke=0.2,shape=21)+
    theme_bw()+
    ylab("allele1&2")+
    ggtitle(ifelse(chrgroupe=="chr1:3",sample," "))+
    theme(axis.title.x = element_blank(),axis.text = element_blank(),
          axis.ticks = element_blank(),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank(),
          strip.background = element_rect(fill="white"),axis.title.y = element_text(size = 8),
          plot.title = element_text(size=10,hjust = 0, vjust = 1))+
    scale_colour_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))+
    scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))
  gridExtra::rbind.gtable(.plt)
}

plot_1kg = mafs %>>%
  group_by(sample,groupe) %>>%
  nest(.key=variations) %>>% #head(5) %>>%
  by_row(.to="plot",plot_by_sample)

ggsave("../test.pdf",gridExtra::marrangeGrob(plot_1kg$plot,nrow = 10,ncol = 1,top = NULL),width = 8,height = 12)




#=======================
# by genes
#=======================
sample_posi = plot_posi %>>%
  mutate(t_=1) %>>%
  mutate(row_num = cumsum(t_)) %>>%
  dplyr::select(gene_symbol,row_num,role)

plot_gposi=maf_list%>>%head(210) %>>%
  filter(allele=="allele2") %>>%
  mutate(t_=1) %>>%
  mutate(posi=cumsum(t_)-1) %>>%
  mutate(groupe= posi %/% 21, gposi = posi %%21) %>>%
  mutate(groupe=paste("groupe",groupe,sep=""))

mafs=maf_list %>>%head(210) %>>%
  mutate(maf=purrr::map(path,~strip_maf(.))) %>>%
  select(-path) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  mutate(gene_width=gene_end - gene_start) %>>%
  filter(start>gene_start,start<gene_end) %>>%
  mutate(freqposi=(start - gene_start)/gene_width) %>>%
  mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi)) %>>%
  dplyr::select(sample,allele,gene_symbol,freqposi,mutype,PolyPhen) %>>%
  left_join(plot_gposi%>>%dplyr::select(-allele),by=c("sample")) %>>%
  mutate(freqposi=freqposi + gposi) %>>%
  dplyr::select(-gposi) %>>%
  mutate(allele=ifelse(allele=="allele1",0,1))

.plot_gposi=plot_gposi %>>%　mutate(focal=1) %>>%left_join(data_frame(row_num=c(1:105))%>>%mutate(focal=1)) %>>%
left_join(sample_posi) %>>%dplyr::select(sample,groupe,gposi,gene_symbol,row_num,role) %>>% filter(gene_symbol=="BRCA1",groupe=="groupe0")%>>%group_by(gene_symbol,groupe)%>>%nest(.key=variations)
plot_by_gene=function(.plot_gposi){
  sample_groupe=.plot_gposi$groupe
  gene=.plot_gposi$gene_symbol
  role=first(.plot_gposi$variations[[1]]$role)
  .maf=mafs %>>% filter(gene_symbol==gene,groupe==sample_groupe) %>>%
    mutate(quoti=(freqposi %/% 0.05)+0.025,surplus= freqposi %% 0.05) %>>%
    group_by(allele,quoti) %>>%
    mutate(rank=dplyr::row_number(surplus)) %>>%
    ungroup() %>>%
    mutate(allele = allele+((rank-1)*0.12+0.04),freqposi=quoti / 20 +0.025) %>>%
    dplyr::select(-quoti,-surplus) 
  .plt=.plot_gposi$variations[[1]] %>>%mutate(groupe=sample_groupe)%>>%
    ggplot()+
    geom_rect(aes(xmin=gposi,xmax=gposi+1,ymin=0,ymax=1),fill="gray")+
    geom_rect(aes(xmin=gposi,xmax=gposi+1,ymin=1,ymax=2),fill="gray50")+
    facet_grid(groupe~.)+
    coord_cartesian(ylim = c(-0.05,2),xlim=c(-0.05,21),expand=F)+
    geom_hline(yintercept = c(1,2),colour="white")+
    geom_vline(xintercept = 0:21,colour="black",size=0.1)+
    geom_point(data=.maf,aes(x=freqposi,y= allele,fill=mutype),size=1,stroke=0.2,shape=21)+
    theme_bw()+
    ylab("allele1&2")+
    xlab("by_sample")+
    ggtitle(ifelse(sample_groupe=="groupe0",paste(gene,role,sep=" : ")," "))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank(),
          strip.background = element_rect(fill="white"),axis.title.y = element_text(size = 8),
          plot.title = element_text(size=20,hjust = 0, vjust = 1))+
    scale_colour_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))+
    scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))
  gridExtra::rbind.gtable(.plt)
}

plot_1kg_bygene = plot_gposi %>>%
  mutate(focal=1) %>>%left_join(data_frame(row_num=c(1:105))%>>%mutate(focal=1)) %>>%
  left_join(sample_posi) %>>%
  dplyr::select(sample,groupe,gposi,gene_symbol,row_num,role) %>>%
  group_by(gene_symbol,row_num,groupe) %>>%
  nest(.key='variations') %>>% dplyr::arrange(row_num,groupe) %>>% #head(5) %>>%
  by_row(.to="plot",plot_by_gene)

ggsave("by_gene_105sample.pdf",gridExtra::marrangeGrob(plot_1kg_bygene$plot,nrow = 10,ncol = 1,top = NULL),width = 8,height = 12)
