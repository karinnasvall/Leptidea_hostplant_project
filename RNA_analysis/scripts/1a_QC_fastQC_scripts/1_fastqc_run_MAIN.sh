#!/bin/bash

###### START UPPMAX RUNS 

echo
echo "Starting Uppmax jobs ..."
echo




##### Paths and folders

#remember initial path
SRCDIR_INI=$(pwd) 

# Path to folder containing input data
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/data/0_RAW_DATA

# Path to results folder
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1_fastqc

# Library list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt


                   






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












##### Start jobs

cd $SRCDIR_INI

for i in `seq 1 1 $N_LIBS`; do                  #loop starting all jobs
   DNAFILE_A=${READS[$i,1]}                  
   DNAFILE_B=${READS[$i,2]}

   GO_SAMPLE=${SAMPLE_LIST[$i]}

   echo $DNAFILE_A                                             #display fastq file names (paired)
   echo $DNAFILE_B
   echo $GO_SAMPLE

   sbatch 1_fastqc_run_query.sh $AA/$DNAFILE_A \
                              $AA/$DNAFILE_B \
                              ${DNAFILE_A} \
                              ${DNAFILE_B} \
                              $RR \
                              ${GO_SAMPLE}


                                                               
   sleep 1                                                     
   echo
    
done


squeue -u $USER                                    #check job status
echo





