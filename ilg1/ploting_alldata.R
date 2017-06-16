#this script is ~/ascat/script/ploting_alldaat.R
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

setwd('/Volumes/areca42TB/tcga/CNA/breast/cel/')
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
topdriver_bed=read_tsv("../../../maf_norm/top_driver105.bed",col_names = c("chr","gene_start","gene_end","ids","score","strand"))%>>%
  tidyr::separate(ids,c("gene_symbol","ensg"),sep=";")
brca_maf=read_tsv("/working/maf/7a07a833-4eab-44c9-bbf6-a64bd51c012e/TCGA.BRCA.muse.7a07a833-4eab-44c9-bbf6-a64bd51c012e.protected.maf.gz",comment = "#")
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


cna=read_tsv('annotate_ascat.tsv.gz') 
ascat_focal=read_tsv('error_of_annotate_ascat.txt') %>>%
  mutate(focal_ascat="no")
maf_focal=read_tsv('../../../maf_norm/breast/depth/list_of_perfect_maf.tsv') %>>%
  left_join(ascat_focal) %>>%
  mutate(focal=ifelse(is.na(focal_ascat),focal,"no")) %>>%dplyr::select(-focal_ascat)
cna = cna %>>%
  left_join(maf_focal) %>>%
  filter(focal=="ok") %>>%
  dplyr::select(-focal)


###========================
# by patient
###==========================
groups=data_frame(group=c(0:4),groupe=c("chr1:3","chr3:8","chr9:13","chr13:17","chr18:X"))
plot_posi = cna%>>%
  filter(str_detect(patient_id,'TCGA-3C-AAAU')) %>>%
  mutate(num=1) %>>%
  mutate(lnum=cumsum(num)-1) %>>%
  mutate(group=lnum %/% 21, gposi=lnum %% 21) %>>%
  left_join(groups)%>>%
  dplyr::select(gene_symbol,groupe,gposi) %>>%
  left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene"))
write_df_watal(plot_posi,"~/ascat/data/plot_posi.tsv")


plot=function(.cna){
  patient=.cna$patient_id
  chrgroupe=.cna$groupe
  norm_maf1=read_tsv(paste('../../../maf_norm/breast/',patient,'.maf',sep=""),comment = "#") %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position) %>>%
    dplyr::select(gene_symbol,chr,start,end,Consequence) %>>%
    classify_consequence() %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
    filter(!is.na(strand)) %>>%
    mutate(gene_width=gene_end - gene_start) %>>%
    filter(start>gene_start,end<gene_end) %>>%
    mutate(freqposi=(start - gene_start)/gene_width) %>>%
    filter(!mutype=="silent") %>>%
    filter(!mutype=="flank") %>>%
    mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi)) %>>%
    dplyr::select(gene_symbol,freqposi,mutype) %>>%
    left_join(plot_posi) %>>%
    filter(groupe==chrgroupe) %>>%
    mutate(freqposi=freqposi + gposi) %>>%
    dplyr::select(-gposi)

  maf1=brca_maf %>>%
    filter(str_detect(Tumor_Sample_Barcode,patient)) %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position) %>>%
    dplyr::select(gene_symbol,chr,start,end,Consequence,t_depth,t_ref_count,t_alt_count) %>>%
    classify_consequence() %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
    filter(!is.na(strand)) %>>%
    mutate(gene_width=gene_end - gene_start) %>>%
    filter(start>gene_start,end<gene_end) %>>%
    mutate(freqposi=(start - gene_start)/gene_width) %>>%
    filter(!mutype=="silent") %>>%
    filter(!mutype=="flank") %>>%
    left_join(.cna$CNA[[1]]%>>%dplyr::select(gene_symbol,nmajor,nminor,purity)%>>%mutate(nall=nmajor+nminor)) %>>%
    mutate(t_ref_count = ifelse((t_ref_count - (1-purity)* t_depth) < 0, 0, t_ref_count - (1-purity)* t_depth)) %>>%
    mutate(alt_allele_count=nall*(t_alt_count/(t_ref_count + t_alt_count))) %>>%
    mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi),
           y_start=ifelse(abs(nmajor - alt_allele_count) >= abs(nminor - alt_allele_count), 0, nminor )) %>>%
    mutate(y_end= y_start + alt_allele_count) %>>%
    dplyr::select(gene_symbol,freqposi,mutype,y_start,y_end) %>>%
    left_join(plot_posi) %>>%
    filter(groupe==chrgroupe)%>>%
    mutate(freqposi=freqposi + gposi) %>>%
    dplyr::select(-gposi)
  
  .plt = .cna$CNA[[1]] %>>%
    ggplot()+
    geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=0,ymax=nminor),fill="gray")+
    geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=nminor,ymax=nminor+nmajor),fill="gray50")+
    facet_grid(groupe~.)+
    coord_cartesian(ylim = c(-0.4,5.2),xlim=c(0,21),expand=F)+
    geom_hline(yintercept = c(1,2),colour="white")+
    geom_vline(xintercept = 0:20,colour="black",size=0.1)+
    annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gene_symbol},
             x=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gposi} + 0.5,
             y=-0.15, size = 1.5,colour="red")+
    annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gene_symbol},
             x=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gposi} + 0.5,
             y=-0.15, size = 1.5,colour="black")+
    geom_segment(data=maf1%>>%filter(groupe==chrgroupe),
                 aes(x=freqposi, y=y_start, xend=freqposi, yend=y_end, colour=mutype))+
    geom_dotplot(data=norm_maf1%>>%filter(groupe==chrgroupe),aes(x=freqposi,fill=mutype,colour=mutype),
                 stackgroups = T, binwidth = 0.05,dotsize=1.5,binpositions = "all")+
    theme_bw()+
    ylab("each allele number")+
    ggtitle(ifelse(chrgroupe=="chr1:3",patient," "))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank(),
          strip.background = element_rect(fill="white"),axis.title.y = element_text(size = 8),
          plot.title = element_text(size=10,hjust = 0, vjust = 1))+
    scale_colour_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))+
    scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))
  gridExtra::rbind.gtable(.plt)
}

pa=cna %>>%
  left_join(plot_posi) %>>%
  mutate(start_rate=start_rate + gposi,end_rate=end_rate + gposi) %>>%
  group_by(patient_id,groupe) %>>%
  nest(.key=CNA) %>>%head(5)%>>%
  by_row(.to='plot',plot)
ggsave("by_patient.pdf",gridExtra::marrangeGrob(pa$plot,nrow = 10,ncol = 1,top = NULL),width = 8,height = 12)

#############test polts#######################
if(1){
  chrgroupe="chr1:3"
  patient="TCGA-3C-AAAU"
  norm_maf1=read_tsv(paste('../../../maf_norm/breast/',patient,'.maf',sep=""),comment = "#") %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position) %>>%
    dplyr::select(gene_symbol,chr,start,end,Consequence) %>>%
    classify_consequence() %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
    filter(!is.na(strand)) %>>%
    mutate(gene_width=gene_end - gene_start) %>>%
    filter(start>gene_start,end<gene_end) %>>%
    mutate(freqposi=(start - gene_start)/gene_width) %>>%
    filter(!mutype=="silent") %>>%
    filter(!mutype=="flank") %>>%
    mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi)) %>>%
    dplyr::select(gene_symbol,freqposi,mutype) %>>%
    left_join(plot_posi) %>>%
    filter(groupe==chrgroupe) %>>%
    mutate(freqposi=freqposi + gposi) %>>%
    dplyr::select(-gposi) %>>%
    mutate(quoti=(freqposi %/% 0.05)+0.025,surplus= freqposi %% 0.05) %>>%
    group_by(quoti) %>>%
    mutate(rank=dplyr::row_number(surplus)) %>>%
    ungroup() %>>%
    mutate(height=((rank-1)*0.12+0.04),freqposi=quoti / 20 +0.025) %>>%
    dplyr::select(-quoti,-surplus)
  
  maf1=brca_maf %>>%
    filter(str_detect(Tumor_Sample_Barcode,patient)) %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position) %>>%
    dplyr::select(gene_symbol,chr,start,end,Consequence,t_depth,t_ref_count,t_alt_count) %>>%
    classify_consequence() %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
    filter(!is.na(strand)) %>>%
    mutate(gene_width=gene_end - gene_start) %>>%
    filter(start>gene_start,end<gene_end) %>>%
    mutate(freqposi=(start - gene_start)/gene_width) %>>%
    filter(!mutype=="silent") %>>%
    filter(!mutype=="flank") %>>%
    left_join(cna%>>%filter(patient_id==patient)%>>%dplyr::select(gene_symbol,nmajor,nminor,purity)%>>%mutate(nall=nmajor+nminor)) %>>%
    mutate(t_ref_count = ifelse((t_ref_count - (1-purity)* t_depth) < 0, 0, t_ref_count - (1-purity)* t_depth)) %>>%
    mutate(alt_allele_count=nall*(t_alt_count/(t_ref_count + t_alt_count))) %>>%
    mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi),
           y_start=ifelse(abs(nmajor - alt_allele_count) >= abs(nminor - alt_allele_count), 0, nminor )) %>>%
    mutate(y_end= y_start + alt_allele_count) %>>%
    dplyr::select(gene_symbol,freqposi,mutype,y_start,y_end) %>>%
    left_join(plot_posi) %>>%
    filter(groupe==chrgroupe)%>>%
    mutate(freqposi=freqposi + gposi) %>>%
    dplyr::select(-gposi)
  cna%>>%filter(patient_id==patient) %>>%
    left_join(plot_posi) %>>%
    mutate(start_rate=start_rate + gposi,end_rate=end_rate + gposi) %>>%
    left_join(driver_genes%>>%dplyr::select(gene,role),by=c("gene_symbol"="gene")) %>>%
    filter(groupe==chrgroupe) %>>%
    ggplot()+
      geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=0,ymax=nminor),fill="gray")+
      geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=nminor,ymax=nminor+nmajor),fill="gray50")+
      facet_grid(groupe~.)+
      coord_cartesian(ylim = c(-0.3,5.2),xlim=c(0,21),expand=F)+
      geom_hline(yintercept = c(1,2),colour="white")+
      geom_vline(xintercept = 0:20,colour="black",size=0.1)+
      annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gene_symbol},
               x=plot_posi%>>%filter(groupe==chrgroupe,role=="TSG")%>>%{.$gposi} + 0.5,
               y=-0.15, size = 1.5,colour="red")+
      annotate("text",label=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gene_symbol},
               x=plot_posi%>>%filter(groupe==chrgroupe,role=="oncogene")%>>%{.$gposi} + 0.5,
               y=-0.15, size = 1.5,colour="black")+
      geom_segment(data=maf1%>>%filter(groupe==chrgroupe),
                   aes(x=freqposi, y=y_start, xend=freqposi, yend=y_end, colour=mutype))+
      geom_point(data=norm_maf1 %>>%filter(groupe==chrgroupe),aes(x=freqposi, y=height, fill=mutype),size=1,stroke=0.2,shape=21)+
#     geom_dotplot(data=norm_maf1%>>%filter(groupe==chrgroupe),aes(x=freqposi),
#                   stackgroups = T, binwidth = 0.05,dotsize=1.6,stackratio = 0.9375,binpositions = "all")+
#     geom_dotplot(data=norm_maf1%>>%filter(groupe==chrgroupe),aes(x=freqposi,colour=mutype,fill=mutype),
#                   stackgroups = T, binwidth = 0.05,dotsize=1.5,binpositions = "all")+
      theme_bw()+
      ylab("each allele number")+
      ggtitle(ifelse(chrgroupe=="chr9:13",patient," "))+
      theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),panel.border = element_blank(),
            strip.background = element_rect(fill="white"),axis.title.y = element_text(size = 8),
            plot.title = element_text(size=20,hjust = 0, vjust = 0.5,angle = 90))+
      scale_colour_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))+
      scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))
}

##==============================
## by gene
##==============================
sample_posi = plot_posi %>>%
  mutate(t_=1) %>>%
  mutate(row_num = cumsum(t_)) %>>%
  dplyr::select(gene_symbol,row_num,role)

maf_focal_bygene=maf_focal %>>%
  filter(focal=="ok") %>>% head(105)

plot_gposi = maf_focal_bygene %>>%
  dplyr::select(-focal) %>>%
  mutate(t_=1) %>>%
  mutate(posi=cumsum(t_)-1) %>>%
  mutate(groupe= posi %/% 21, gposi = posi %%21) %>>%
  mutate(groupe=paste("groupe",groupe,sep=""))

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

norm_mafs = maf_focal_bygene%>>%
  mutate(filename=paste('../../../maf_norm/breast/',patient_id,'.maf',sep="")) %>>%
  mutate(purrr::map(filename,~strip_maf(.))) %>>%
  unnest() %>>%
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position) %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  mutate(gene_width=gene_end - gene_start) %>>%
  filter(start>gene_start,start<gene_end) %>>%
  mutate(freqposi=(start - gene_start)/gene_width) %>>%
  mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi)) %>>%
  dplyr::select(patient_id,gene_symbol,freqposi,mutype) %>>%
  left_join(plot_gposi) %>>%
  mutate(freqposi=freqposi + gposi) %>>%
  dplyr::select(-gposi) %>>%
  mutate(quoti=(freqposi %/% 0.05)+0.025,surplus= freqposi %% 0.05) %>>%
  group_by(gene_symbol,patient_id,quoti) %>>%
  mutate(rank=dplyr::row_number(surplus)) %>>%
  ungroup() %>>%
  mutate(allele =((rank-1)*0.2+0.04),freqposi=quoti / 20 +0.025) %>>%
  dplyr::select(-quoti,-surplus)

brca_maf_bygene = brca_maf %>>%
  mutate(patient_id=str_extract(Tumor_Sample_Barcode,"TCGA-[^-]*-[^-]*"))%>>%
  left_join(maf_focal_bygene) %>>% filter(!is.na(focal)) %>>% 
  dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position) %>>%
  dplyr::select(gene_symbol,chr,start,end,Consequence,t_depth,t_ref_count,t_alt_count,patient_id,PolyPhen) %>>%
  classify_consequence() %>>%
  left_join(topdriver_bed%>>%dplyr::select(gene_symbol,gene_start,gene_end,strand)) %>>%
  filter(!is.na(strand)) %>>%
  mutate(gene_width=gene_end - gene_start) %>>%
  filter(start>gene_start,end<gene_end) %>>%
  mutate(freqposi=(start - gene_start)/gene_width) %>>%
  filter(!mutype=="silent") %>>%
  filter(!mutype=="flank") %>>%
  left_join(cna%>>%dplyr::select(gene_symbol,nmajor,nminor,purity,patient_id)%>>%mutate(nall=nmajor+nminor),
            by=c("gene_symbol","patient_id")) %>>%
  mutate(t_ref_count = ifelse((t_ref_count - (1-purity)* t_depth) < 0, 0, t_ref_count - (1-purity)* t_depth)) %>>%
  mutate(alt_allele_count=nall*(t_alt_count/(t_ref_count + t_alt_count))) %>>%
  mutate(freqposi=ifelse(strand=="-",1-freqposi,freqposi),
         y_start=ifelse(abs(nmajor - alt_allele_count) >= abs(nminor - alt_allele_count), 0, nminor )) %>>%
  mutate(y_end= y_start + alt_allele_count) %>>%
  dplyr::select(gene_symbol,patient_id,freqposi,mutype,y_start,y_end,PolyPhen) %>>%
  left_join(plot_gposi) %>>%
  mutate(freqposi=freqposi + gposi) %>>%
  dplyr::select(-gposi)

cna_bygene = cna %>>%
  left_join(maf_focal_bygene) %>>%filter(!is.na(focal)) %>>%
  left_join(plot_gposi) %>>%
  mutate(start_rate=start_rate + gposi,end_rate=end_rate + gposi) %>>%
  dplyr::select(gene_symbol,groupe,patient_id,nmajor,nminor,start_rate,end_rate)

.cna=cna_bygene %>>%filter(gene_symbol=="PHF6",groupe=="groupe0")%>>%mutate(role="TSG") %>>%group_by(gene_symbol,groupe) %>>%nest(.key=CNA)
plot_by_gene=function(.cna){
  gene=.cna$gene_symbol
  patient_groupe=.cna$groupe
  role=first(.cna$CNA[[1]]$role)
#  print(paste(gene,patient_groupe,sep=" "))
  .norm_maf=norm_mafs %>>%filter(gene_symbol==gene,groupe==patient_groupe)
  .maf=brca_maf_bygene %>>%filter(gene_symbol==gene,groupe==patient_groupe)
  .plt=.cna$CNA[[1]] %>>%
    ggplot()+
    geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=0,ymax=nminor),fill="gray")+
    geom_rect(aes(xmin=start_rate,xmax=end_rate,ymin=nminor,ymax=nminor+nmajor),fill="gray50")+
    coord_cartesian(ylim = c(-0.05,5.2),xlim=c(-0.05,21),expand=F)+
    geom_hline(yintercept = c(1,2),colour="white")+
    geom_vline(xintercept = 0:21,colour="black",size=0.1)+
    geom_segment(data=.maf,aes(x=freqposi, y=y_start, xend=freqposi, yend=y_end, colour=mutype))+
    geom_point(data=.norm_maf,aes(x=freqposi, y=allele, fill=mutype),size=0.8,stroke=0.1,shape=21)+
#    facet_grid(groupe~.)+
    theme_bw()+
    ylab("each allele number")+
    ggtitle(ifelse(patient_groupe=="groupe0",paste(gene,role,sep=" : ")," "))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank(),
          strip.background = element_rect(fill="white"),axis.title.y = element_text(size = 8),
          plot.title = element_text(size=20,hjust = 0, vjust = 0.5))+
    scale_colour_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))+
    scale_fill_manual(values=c( missense="forestgreen",truncating="black",splice="blue"))
  gridExtra::rbind.gtable(.plt)
}

plot_bygene = cna_bygene %>>%
  left_join(sample_posi) %>>%
  group_by(gene_symbol,row_num,groupe)%>>%nest(.key=CNA) %>>%dplyr::arrange(row_num,groupe) %>>%
  by_row(.to="plot",~plot_by_gene(.))

ggsave("plot_bygene.pdf",gridExtra::marrangeGrob(plot_bygene$plot,nrow = 10,ncol = 1,top = NULL),width = 8,height = 12)
