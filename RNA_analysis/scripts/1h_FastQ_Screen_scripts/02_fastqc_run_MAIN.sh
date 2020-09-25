#!/bin/bash

###### START UPPMAX RUNS

echo
echo "Starting Uppmax jobs ..."
echo




##### Paths and folders

# Path to folder containing input data
AA=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1h_FastQ_Screen

# Path to results folder
RR=$AA

# Library list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt





#remember initial path
SRCDIR_INI=$(pwd)






######## read sample list

i=1

while read -r LINE
do
   SAMPLE_LIST[$i]=$LINE
   let "i+=1"
done < ${SPL}


#Number of libraries
N_LIBS=${#SAMPLE_LIST[@]}





######## file names


declare -A READS

for i in `seq 1 1 $N_LIBS`; do
   READS[$i,1]=${SAMPLE_LIST[$i]}_R1_001_FASTQScreen_filtered.fq.gz
   READS[$i,2]=${SAMPLE_LIST[$i]}_R2_001_FASTQScreen_filtered.fq.gz
done








#cd $PATH_MAIN/                                     #cd to project folder
#OUT_FOLDER=$SRCDIR_INI                                 #output folder
#if [ -d "$OUT_FOLDER" ]                            #check whether output folder already exists
#   then
#      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
#   else
#      mkdir $OUT_FOLDER                            #create output folder
#fi




##### Start jobs

cd $SRCDIR_INI

for i in `seq 1 1 $N_LIBS`; do                  #loop starting all jobs
   DNAFILE_A=${READS[$i,1]}
   DNAFILE_B=${READS[$i,2]}

   GO_SAMPLE=${SAMPLE_LIST[$i]}

   echo $DNAFILE_A                                             #display fastq file names (paired)
   echo $DNAFILE_B
   echo $GO_SAMPLE

   sbatch 02_fastqc_run_query.sh $AA/$DNAFILE_A \
                                 $AA/$DNAFILE_B \
                                 ${DNAFILE_A} \
                                 ${DNAFILE_B} \
                                 $RR \
                                 ${GO_SAMPLE}


                                                               #note: "%.fastq.gz" removes the ".fastq.gz" suffix
    sleep 1                                                     #pauses for 1 sec
    echo
done


squeue -u $USER                                    #check job status
echo




