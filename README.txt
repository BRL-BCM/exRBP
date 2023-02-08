#Counting genes will use the hg19 gtf file and information in
#Biotype calculations for the waffle plot figures also use this information
/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/Cell_line_all/Intersect_with_RefSeq

#Begin by getting the bedgraphs for the studies of interest
#No changes need to be made to this script other than updating the totalN to the studies you want
launch_study_loop.R
create_loci_by_study_matrixes_loop.pbs
create_study_matrix_sum.sh # This now collapses each of the 150 RBPs down to 1 column    
#If study is too large or you need to only process some of the biosamples, use these as examples. It tests each biosample against a file list before processing it and adding it to the output file
qsub -N "LLAUR_test" -v study="EXR-LLAUR1SFszrF-AN" create_study_matrix_llaur.pbs
create_study_matrix_llaur.pbs
create_study_matrix_sum_llaur.sh

#Then if you are not using all samples within that study, filter them to the samples you want
#make sure to accurately update the output file names as you will use that in the next step to search for all files related to the analysis you are doing
launch_filter_loop.R
filter_full_study.pbs #If using collapsed version, use - and if using original version use . for sample names (EXR-/EXR.)
#Check this file for the list which samples to filter for for each biofluid
biofluid_samples.txt
#If you are using all samples, simply copy the files so that they all have the same suffix. You will need it so you can search for all files related to the analysis you are trying to do in the next step.
filter_full_study_cellLines.pbs

#Then combine all the studies you are interested in into one file
qsub -v type=cellLine,status=all -N combined_filter_cellLine combined_and_filter.pbs
qsub -v type=serum,status=healthy -N combine_serum combined_and_filter.pbs
#combined_and_filter_pht.pbs
#You will need to update the file you are searching for and make sure to update the output file names


####THIS STEP HAS BEEN REMOVED WITH THE ADDITION OF COLLAPSING FILES IN STEP 1
#Finally, combine all the individual columns (one for each RBP intersect) into one column for each file
#collapse_columns.pbs
#Make sure to update to the correct file name. You will also need to give it the name of the first biosample in the file so be sure to open the study file and look at what biosample is first.


####To find bases covered in exRNA
#Collpase each study using unionbg
launch_all_exRNA.R
all_exRNA_coverage.pbs
#Collapse down to number of samples per region that have a given coverage
launch_all_exRNA_collapse_large.R #In here, set the coverage, coverage=...
all_exRNA_collapse_largeStudy3.pbs
#This one get the coverage for any site that has reads
exRNA_coverage.sh
#Combine across studies using unionbg
```
module load BEDTools/2.17
sed -i $'s/ /\t/g' *_coverage_1.bedgraph
sed -i '1d' *_coverage_1.bedgraph
cov=1
locations=$(find /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/study_exRNA_coverage -type f -name '*_coverage_1.bedgraph' | sort |sed -e 's/\n/\t/g')
names=$(ls /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/study_exRNA_coverage/*_coverage_${cov}.bedgraph | sort | awk -F_coverage_${cov}.bedgraph '{print $1}' | cut -c 72-)
bedtools unionbedg -i $locations -header -names $names > /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/study_exRNA_coverage/studies_intersected_${cov}.bed
#Column  1-3 are region and column 4 is subtracting end from stop to get length of each region
awk -v OFS="\t" '{print $1,$2,$3}' studies_intersected_${cov}.bed | awk -v OFS="\t" '{print $1,$2,$3,$3-$2}'>  studies_intersected_${cov}_regions.bed
awk -F'\t' '{sum+=$4;} END{print sum;}' studies_intersected_${cov}_regions.bed #sums column4 (which is the length of each region) to get coverage of genome
gzip studies_intersected_${cov}.bed
gzip studies_intersected_${cov}_regions.bed
```

###To then find those that overlap with the eCLIP
#eCLIP_coverage.txt is the first 3 columns of the intersect150.bed file for rows that have at least 1 hit
awk -F"\t" '$3>$2' eCLIP_coverage.txt | sort -k1,1 -k2,2n > eCLIP.bed 
cov=2
for FILE_i in *_${cov}.bedgraph
 do
  PREFIX=${FILE_i::${#FILE_i}-20}
  echo "${PREFIX}"
  sort -k1,1 -k2,2n $FILE_i > tmp
  #For each study map how many samples had reads onto eCLIP locations
  bedtools map -a eCLIP.bed -b tmp -c 4 | sed -e 's/\./0/g' > ${PREFIX}_${cov}_eCLIP.bedgraph
#  wc -l ${PREFIX}_${cov}_eCLIP.bedgraph
  rm tmp
 done

#Combine all studies file then keep the first 3 location columns and the read columns for every study
paste *_${cov}_eCLIP.bedgraph | awk -v OFS="\t" '{print $1, $2, $3, $4, $8,$12, $16, $20, $24, $28, $32, $36, $40, $44, $48, $52, $56, $60, $64, $68,$72, $76, $80, $84, $88, $92, $96, $100, $104, $108, $112, $116, $120, $124,$128, $132, $136, $140, $144, $148, $152, $156, $160, $164, $168, $172, $176,$180, $184, $188, $192, $196, $200, $204, $208, $212, $216, $220, $224, $228,$232, $236, $240, $244, $248, $252, $256, $260, $264, $268, $272}' > studies_intersected_${cov}_eCLIP.bed
#Count number of studies where a sample has hits keep rows over 1
awk -v COV=1 -v OFS="\t" '{check=0;for(i=4;i<=NF;i++)if($i>=COV)check+=1;print $1,$2,$3,check}' studies_intersected_${cov}_eCLIP.bed > tmp
awk -F"\t" 'NR==1{print;next}$4>0' tmp > studies_intersected_${cov}_eCLIP_filtered.bed
#Find length of each region where a sample hit and print sum
awk -v OFS="\t" '{print $1,$2,$3}' studies_intersected_${cov}_eCLIP_filtered.bed | awk -v OFS="\t" '{print $1,$2,$3,$3-$2}'>  studies_intersected_${cov}_eCLIP_regions.bed
awk -F'\t' '{sum+=$4;} END{print sum;}' studies_intersected_${cov}_eCLIP_regions.bed
rm tmp



