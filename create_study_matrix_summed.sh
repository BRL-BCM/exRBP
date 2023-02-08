#!/usr/bin/bash

#STUDY_LOC="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EXR-KJENS10lPClY-AN"
STUDY_LOC="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EXR-DGALA1QDi9GG-AN"
cd $STUDY_LOC

STUDY=$(basename $STUDY_LOC)
if [ ! -d /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/${STUDY} ]; then
  mkdir -p /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/${STUDY};
fi

for biosample in ${STUDY_LOC}/RBP_intersects/biosample_RBP/*.xz
do
 #If study hasn't finished processing so the RBP_intersects folder or the biosample_RBP doesn't exist skip to next study
 [ -f "$biosample" ] || break
 cd ${STUDY_LOC}/RBP_intersects/biosample_RBP/
 #The name comes with the path attached so create the prefix which is just the biosample name
 PREFIX_i=$(basename $biosample)
 echo ${PREFIX_i}
 SHORT=${PREFIX_i%.*}
 
 #If biosample file doesn't exist, create it
 if [ ! -f /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/${STUDY}/${SHORT} ]; then
  #unzip
  xz -k -d $biosample
  
  mv $SHORT /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/${STUDY}/${SHORT}
 fi
 
done

cd /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/${STUDY}
files=$(ls *.bed | head -1)

for bedfile in *.bed
do
 test=$(basename $bedfile)
 #echo $test
 SAMPLE=$(basename $bedfile | sed 's/-BS.*/-BS_/')
 #Get rid of unwanted columns
 #awk '{print $4, $5, $6, $7}' $bedfile | sed 's/ /\t/g' > sample_total_1_temp.txt #short for testing
 #Select only columns where RBP has been in 1 cell line or the merged file if it has been in 2. No need to take individual cell lines if merged file is present
 awk '{print $4, $5, $6, $7, $10, $13, $14, $15, $18, $19, $20, $23, $24, $25, $26, $29, $30, $31, $32, $35, $36, $37, $40, $41, $42, $43, $46, $49, $52, $53, $56, $59, $60, $61, $62, $63, $64, $67, $70, $73, $74, $75, $78, $79, $82, $83, $86, $87, $88, $89, $90, $91, $94, $97, $100, $103, $106, $109, $112, $115, $118, $121, $124, $125, $126, $129, $130, $133, $136, $139, $142, $145, $148, $149, $150, $153, $154, $155, $156, $157, $160, $161, $162, $163, $164, $165, $168, $169, $170, $171, $172, $173, $174, $177, $180, $181, $182, $183, $186, $189, $192, $195, $196, $197, $200, $201, $204, $205, $208, $209, $210, $211, $214, $215, $216, $219, $222, $225, $228, $231, $232, $235, $236, $237, $238, $241, $244, $245, $248, $251, $252, $255, $258, $261, $264, $267, $270, $273, $274, $275, $278, $279, $280, $283, $286, $289, $290, $293, $294, $295, $298}' $bedfile | sed 's/ /\t/g' > sample_total_1_temp.txt
 awk '{c=0;for(i=1;i<=NF;++i){c+=$i};print c}' sample_total_1_temp.txt > sample_total_1.txt
 #concatenate
 if [ $test == $files ]; then
  #echo "equal"
  mv sample_total_1.txt study_total_1.txt
 else
  #echo "Not"
  paste study_total_1.txt sample_total_1.txt > study_total_1_temp.txt
  mv study_total_1_temp.txt study_total_1.txt
 fi
 echo $SAMPLE > sample_finished.txt
done

#add column names to study
NAMES_USE=$(ls *.bed | sed 's/_EXR-.*.//g')
echo $NAMES_USE | cat - study_total_1.txt > study_total_1_temp.txt && mv study_total_1_temp.txt study_total_1.txt
#concatenate columns 1-3
awk '{print $1, $2, $3}' $bedfile | paste -d "\t" - study_total_1.txt | sed 's/ /\t/g' > ${STUDY}_1.txt
#filter
awk 'NR==1 {print } NR>1 { s=0 ; for(i=4;i<=NF;i++) s+=$i ; if (s) print ;}' ${STUDY}_1.txt > ${STUDY}_filtered1.txt

tar -czvf ${STUDY}_summed.tar.gz ${STUDY}_1.txt
tar -czvf ${STUDY}_summed_filtered.tar.gz ${STUDY}_filtered1.txt

rm ${STUDY}_1.txt
rm ${STUDY}_filtered1.txt
rm *.bed
rm study_total_1.txt
rm sample_total_1.txt
rm sample_total_1_temp.txt
