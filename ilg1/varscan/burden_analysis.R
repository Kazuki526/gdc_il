######！！！！！！！！注意！！！！！！！！########
# prepare_tbl.Rを実行してからこのスクリプトを開始する！！！

#prepare 1000genomes UK10K EXAC(nonTCGA) data sets
.colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2',
              'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position",'Transcript_ID')
.cols = .colnames %>>%
{setNames(c('c','c','d','d','c','c', 'c','c','c','c','c','c'), .)} %>>%
{do.call(readr::cols_only, as.list(.))}
strip_maf = function(infile) {
  read_tsv(infile, comment='#', col_types=.cols) %>>%
    classify_consequence() %>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
    left_join(topdriver_bed%>>%dplyr::select(gene_symbol,strand)) %>>%filter(!is.na(strand))
}
###1kg###

sample_num_1kg=5008
sample_num_X_1kg=sample_num_1kg - 1233
vcf_1kg_ac0=read_tsv("/working/1000genomes/maf/extract/ac0_variation.tsv")
maf_1kg_all=read_tsv("/working/1000genomes/maf/extract/all.maf",col_types = "cccciiccccc")

tally_1kg_all= maf_1kg_all %>>%
  filter(!str_detect(Consequence,"intron_variant")) %>>%
  filter(!mutype=="flank") %>>%
  filter(!str_detect(Consequence, "frame")) %>>%
  group_by(gene_symbol,ref,alt,chr,start,end,mutype,PolyPhen,Consequence) %>>%
  summarise(ac_1kg=n()) %>>% ungroup() %>>%
  left_join(maf_1kg_all %>>% group_by(sample,chr,start,alt) %>>%mutate(variant_num=n()) %>>%
              filter(variant_num==2) %>>%group_by(chr,start,end,alt) %>>%summarise(homo_1kg=n())) %>>%
  mutate(an_1kg = ifelse(chr=="X",sample_num_X_1kg,sample_num_1kg))%>>%
  mutate(chr=paste0("chr",chr),homo_1kg=ifelse(is.na(homo_1kg),0,homo_1kg)) #%>>%
  #full_join(vcf_1kg_ac0)

rm(vcf_1kg_ac0,maf_1kg_all,sample_num_1kg,sample_num_X_1kg)
####UK 10K ###
vcf_10k=read_tsv("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle_likevcf.tsv.gz",col_types = "cdccdddd")
maf_10k=strip_maf("/Volumes/areca42TB/ega/file/all_sample_topdriver_region_indelhandle.maf")
uk10k=maf_10k %>>%
  left_join(vcf_10k%>>%mutate(start=ifelse(ref=="-",start -1,start))) %>>%
  mutate(chr=paste0("chr",chr)) %>>%
  dplyr::rename(ac_uk10k=uk_ac,an_uk10k=uk_an)%>>%
  filter(an_uk10k >3000)
rm(vcf_10k,maf_10k)

### EXAC nonTCGA ###
vcf_exac=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf.tsv.gz",
                  col_types = "cdccdd") %>>%dplyr::rename(start=posi) 
exac_nonindel=strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver.maf") %>>%
  left_join(vcf_exac)
vcf_exac_indel=read_tsv("/Volumes/areca42TB/exac/file/exac_nontcga_liftovered_checked_likevcf_indel.tsv.gz",
                        col_types = "cdccdd") %>>%dplyr::rename(start=posi)
exac_indel = strip_maf("/Volumes/areca42TB/exac/file/exac_nontcga_topdriver_indel.maf") %>>%
  left_join(vcf_exac_indel)
vcf_exac =rbind(vcf_exac,vcf_exac_indel) %>>%mutate(chr=paste0("chr",chr))
exac = rbind(exac_nonindel,exac_indel) %>>%mutate(chr=paste0("chr",chr))
rm(vcf_exac_indel,exac_nonindel,exac_indel)

######################################################################################################
#####################################  classify by ExAC  #############################################
######################################################################################################
classed_site=tally_norm_maf %>>%
  left_join(vcf_exac)%>>%
  left_join(coverage_all ) %>>%
  mutate(AF = ifelse(is.na(ac_exac),0,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))

classed_1kg = tally_1kg_all %>>%
  left_join(vcf_exac) %>>%
  mutate(AF = ifelse(is.na(ac_exac),0,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))
  
classed_uk10k = uk10k %>>%
  left_join(vcf_exac) %>>%
  mutate(AF = ifelse(is.na(ac_exac),0,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))

classed_exac=exac %>>%
  mutate(AF = ifelse(is.na(ac_exac),0,ac_exac/an_exac)) %>>%
  mutate(MAF = ifelse(AF > 0.5, 1 - AF, AF))

########################################################################################################
##################################  remove duplicate region in TCGA  ###################################
########################################################################################################
#### prepare_tbl のerror_region,remove_duplicateを使う
####　あとintronやflankを取り除く
classed_site = classed_site%>>%remove_duplicate() %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant")

classed_1kg = classed_1kg%>>%remove_duplicate() %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") #%>>%
  #filter(!(ac_1kg>130 & homo_1kg <10))

classed_uk10k = classed_uk10k%>>%remove_duplicate() %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant,non_coding_transcript_variant") #%>>%
  #filter(!(uk_hetero >=400 & uk_althomo <= 10))

classed_exac = classed_exac%>>%remove_duplicate() %>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant") %>>%
  filter(Consequence!="intron_variant , non_coding_transcript_variant") 

### transcript
transcript_id=norm_maf %>>%count(gene_symbol,Transcript_ID)%>>%
  dplyr::select(-n)%>>%mutate(focal="ok")%>>%dplyr::select(-focal)


##############################################################################################################
#3つのデータベースで値が大きく違うところ抽出！
if(0){
classed_site %>>%
  mutate(afc = ac_cancer/an_cancer*100)%>>%
  mutate(class_c =ifelse(afc<1,"a",ifelse(afc<5,"b","c")) )%>>%dplyr::select(-AF,-MAF)%>>%
  full_join(classed_1kg %>>%mutate(af1k = ac_1kg/an_1kg*100)%>>%
              mutate(class_1k = ifelse(af1k<1,"a",ifelse(af1k<5,"b","c")))%>>%dplyr::select(-AF,-MAF))%>>%
  full_join(classed_uk10k %>>%mutate(afuk = uk_ac/uk_an*100)%>>%
              mutate(class_uk = ifelse(afuk<1,"a",ifelse(afuk<5,"b","c")))%>>%dplyr::select(-AF,-MAF))%>>%
  full_join(classed_exac %>>%mutate(AF=round(AF*1000000000)/10000000)%>>%
              mutate(class_ex = ifelse(AF<1,"a",ifelse(AF<5,"b","c"))))%>>%
  mutate(class_c=ifelse(is.na(class_c),"a",class_c),class_1k=ifelse(is.na(class_1k),"a",class_1k),
         class_uk=ifelse(is.na(class_uk),"a",class_uk),class_ex=ifelse(is.na(class_ex),"a",class_ex))%>>%
  filter(class_c=="c"|class_uk=="c"|class_uk=="c"|class_ex=="c"|
           !(class_c==class_1k & class_1k==class_uk | class_uk==class_ex))%>>%
  dplyr::select(gene_symbol,chr,start,end,ref,alt,mutype,an_cancer,an_cancer,homo_cancer,afc,an_1kg,homo_1kg,af1k,
                uk_an,uk_althomo,afuk,an_exac,AF,class_c,class_1k,class_uk,class_ex) %>>%
  filter(chr!="chrX",mutype!="flank",mutype!="splice_region")%>>%
  arrange(afc)%>>%write_df("~/Dropbox/install/fourdatabase.tsv")
}
#uk10kのan<3000とduplicate

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#ここから解析
#kokokarakaiseki
#here start analysis

write_df(maf_1kg_all,"~/Dropbox/install/tvz/temporary/maf_1kg_all.tsv.gz")
write_df(norm_maf,"~/Dropbox/install/tvz/temporary/norm_maf.gz")
write_df(tally_norm_maf,"~/Dropbox/install/tvz/temporary/tally_norm_maf.gz")
write_df(tally_1kg_all,"~/Dropbox/install/tvz/temporary/tally_1kg_all.gz")
write_df(uk10k,"~/Dropbox/install/tvz/temporary/uk10k.gz")
write_df(exac,"~/Dropbox/install/tvz/temporary/exac.gz")
write_df(vcf_exac,"~/Dropbox/install/tvz/temporary/vcf_exac.gz")

#####################################################################################################################
### Charlesら(NatCom2015の解析の確認)
TFT_by_row = function(.tbl){
  TFT = function(case_f,control_f,case,control){
    case_m = round(case * case_f)
    control_m = round(control * control_f)
    prop.test(c(case_m,control_m),c(case,control),alternative = "greater")$p.value
  }
  n_can=12704
  n_1kg=5008
  n_uk=4090
  n_exac=106210
  .tbl %>>%
    mutate(to_1kg = ifelse((af_can==0 & af_1kg ==0),1, TFT(af_can,af_1kg ,n_can,n_1kg )),
           to_uk  = ifelse((af_can==0 & af_uk  ==0),1, TFT(af_can,af_uk  ,n_can,n_uk  )),
           to_exac= ifelse((af_can==0 & af_exac==0),1, TFT(af_can,af_exac,n_can,n_exac)))%>>%
    dplyr::select(to_1kg,to_uk,to_exac)
}

total_frequency_truncate = tally_norm_maf %>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_cancer)%>>%
  left_join(varscan_error) %>>% mutate(an_error = ifelse(is.na(an_error),0,an_error)) %>>%
  mutate(start = ifelse(alt == "-",start -1,start)) %>>%
  left_join(coverage_all ) %>>%
  filter(!is.na(an_cancer)) %>>%
  mutate(an_cancer =an_cancer - an_error) %>>%
  dplyr::select(-an_error) %>>%
  full_join(tally_1kg_all%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_1kg,an_1kg))%>>%
  full_join(uk10k%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_uk10k,an_uk10k))%>>%
  full_join(exac%>>%dplyr::select(gene_symbol,chr,start,end,ref,alt,Consequence,mutype,ac_exac,an_exac))%>>%
  filter((mutype=="truncating"|mutype=="splice") & chr!="chrX")%>>%
  remove_duplicate()%>>%
  mutate(af_can  = ifelse(is.na(ac_cancer),0,ac_cancer/an_cancer),
         af_1kg  = ifelse(is.na(ac_1kg   ),0,ac_1kg   /an_1kg   ),
         af_uk   = ifelse(is.na(ac_uk10k ),0,ac_uk10k /an_uk10k ),
         af_exac = ifelse(is.na(ac_exac  ),0,ac_exac  /an_exac  )) %>>%
  filter(af_exac<0.0005)%>>%
  group_by(gene_symbol)%>>%
  summarise(af_can=sum(af_can),af_1kg=sum(af_1kg),af_uk=sum(af_uk),af_exac=sum(af_exac))%>>%
  nest(-gene_symbol) %>>%
  mutate(tft=purrr::map(data,~TFT_by_row(.)))%>>%
  unnest()


