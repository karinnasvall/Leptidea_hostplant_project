#!/bin/bash

###### START UPPMAX RUNS (CALLS trimgalore_run_query.sh BASH FILE FOR EACH SET OF PAIRED-END READS)

echo
echo "Starting Uppmax jobs ..."
echo



##### Paths and folders

# Path to folder containing input data
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/data/0_RAW_DATA

# Path to output folder
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1b_trimgalore

# Library list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt


#Number of libraries
#N_LIBS=${#SAMPLE_LIST[@]}




#remember initial path
SRCDIR_INI=$(pwd)





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


declare -A READS

for i in `seq 1 1 $N_LIBS`; do
   READS[$i,1]=${SAMPLE_LIST[$i]}_L002_R1_001.fastq.gz
   READS[$i,2]=${SAMPLE_LIST[$i]}_L002_R2_001.fastq.gz
done








# prepare output folder

#if [ -d "$OUT_FOLDER" ]                            #check whether output folder already exists
#   then
#      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
#   else
#      mkdir $OUT_FOLDER                            #create output folder
#fi



# start jobs

cd $SRCDIR_INI

for i in `seq 1 1 $N_LIBS`; do
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   GO_SAMPLE=${SAMPLE_LIST[$i]}

   echo $RNAFILE_A                                          #display fastq file names (paired)
   echo $RNAFILE_B
   echo $GO_SAMPLE

   sbatch 01_trimgalore_query.sh $AA/$RNAFILE_A \
                                     $AA/$RNAFILE_B \
                                     ${RNAFILE_A} \
                                     ${RNAFILE_B} \
                                     $RR \
                                     ${GO_SAMPLE}

   sleep 1                                                  #pauses for 1 sec
   echo
done


squeue -u $USER                                    #check job status
echo


