#!/bin/bash


###### START UPPMAX RUNS

#Script used to glean data from summary file produced by FastQC (goes through all zipped files in folder containing fastQC results)

echo
echo "Starting Uppmax job ..."
echo

##### load modules
module load bioinfo-tools
module load python3/3.6.0


#### #paths


#remember current path
SRCDIR=$(pwd)

#input folder
AA=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1g_SortMeRNA

#input file (this file is created from the html output of MultiQC)
FF=multiqc_table.txt

# Sample list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt

#results folder
RR=$SRCDIR

#output file
OUT_file=MRL_content_OUT1_post_SortMeRNA.txt







######## read sample list

i=1

while read -r LINE
do
   SAMPLE_LIST[$i]=$LINE
   let "i+=1"
done < ${SPL}



### Number of libraries
N_LIBS=${#SAMPLE_LIST[@]}






######## file names


declare -A ZFILES

for i in `seq 1 1 $N_LIBS`; do
   ZFILES[$i,1]=${SAMPLE_LIST[$i]}_R1_001_SORTMERNA_CLEAN
   ZFILES[$i,2]=${SAMPLE_LIST[$i]}_R2_001_SORTMERNA_CLEAN
done







##### Get data for each library


for i in `seq 1 1 $N_LIBS`; do                  #loop starting all jobs

   # library info
   zfile1=${ZFILES[$i,1]}
   zfile2=${ZFILES[$i,2]}
   echo
   echo $zfile1
   echo $zfile2


   # glean MRL content
   tns_AUX1="$(grep $zfile1 $AA/$FF | cut -f 5 | cut -f 1 -d ' ')"
   tns_AUX2="$(grep $zfile2 $AA/$FF | cut -f 5 | cut -f 1 -d ' ')"

   # store results
   RES_R1[i]="$zfile1"
   RES_R1[i]="$zfile1"
   RES_tns1[i]="$tns_AUX1"
   RES_tns2[i]="$tns_AUX2"


done





##### save to file



# put results in the right order and save to file
for i in `seq 1 1 $N_LIBS`; do
   if [[ "${i}" = "1" ]]
   then
       echo ${ZFILES[$i,1]} ${RES_tns1[$i]} > $RR/$OUT_file
       echo ${ZFILES[$i,2]} ${RES_tns2[$i]} >> $RR/$OUT_file
   else
       echo ${ZFILES[$i,1]} ${RES_tns1[$i]} >> $RR/$OUT_file
       echo ${ZFILES[$i,2]} ${RES_tns2[$i]} >> $RR/$OUT_file
   fi
done




#runtime:




