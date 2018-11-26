#first of all (mmr.tsv was make by hand)
#perl gene_list2bed_json.pl
#perl gene_list2bed.pl 2
echo "maked bed bed_json"

cd /Volumes/areca42TB2/gdc/mmr

for PROJECT in "brca" "crc" "gbm" "hnsc" "kcc" "lgg" "luad" "lusc" "ov" "prad" "thca" "ucec"
do
perl ~/git/gdc_il/ilg1/varscan/mmr/download.pl $PROJECT
perl ~/git/gdc_il/ilg1/varscan/mmr/varscan.pl $PROJECT
done

perl ~/git/gdc_il/ilg1/varscan/mmr/make_depth_file.pl
perl ~/git/gdc_il/ilg1/varscan/mmr/mk_by_patient_coverage_depth.pl
