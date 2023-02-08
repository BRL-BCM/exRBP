for study in /mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/study_exRNA_coverage/*.bed.tgz
do
# echo ${study}
 tar -xzvf $study
 file_name=${study::${#study}-4}
 PREFIX_in=${study::${#study}-17}
 OUTFILE2="${PREFIX_in}coverage_1.bedgraph"
# echo ${file_name}
# echo ${OUTFILE2}
 awk '{print $1, $2, $3}' $file_name | awk '{print $0, "1"}' > ${OUTFILE2}
 rm $file_name
done

