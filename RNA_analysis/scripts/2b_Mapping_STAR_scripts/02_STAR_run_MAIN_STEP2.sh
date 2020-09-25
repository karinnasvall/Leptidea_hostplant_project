#!/bin/bash -l

###### START UPPMAX RUNS

# Performs read mapping with STAR and transcript quantification with StringTie

echo
echo "Starting Uppmax jobs ..."
echo






##### Paths and folders

# Path to folder containing input data
AA=//proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1g_SortMeRNA

# Path to output folder
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2b_Mapping_STAR_STEP2

# Path to results fron STEP 1
ST1=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2b_Mapping_STAR_STEP1

# Sample list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt

# Reference genome (folder)
REFG=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/LS_genome/star_index/





#remember initial path
SRCDIR_INI=$(pwd)






######## read sample list

i=1

while read -r LINE
do
   SAMPLE_LIST[$i]=$LINE
   let "i+=1"
done < ${SPL}



### Number of samples
N_SAMPLES=${#SAMPLE_LIST[@]}














######## file names


declare -A READS

for i in `seq 1 1 $N_SAMPLES`; do
   READS[$i,1]=${SAMPLE_LIST[$i]}_R1_001_SORTMERNA_CLEAN.fastq.gz
   READS[$i,2]=${SAMPLE_LIST[$i]}_R2_001_SORTMERNA_CLEAN.fastq.gz
done










#start jobs

for i in `seq 1 1 $N_SAMPLES`; do

   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   GO_SAMPLE=${SAMPLE_LIST[$i]}

   echo $RNAFILE_A
   echo $RNAFILE_B
   echo $GO_SAMPLE
   echo ${REFG}

   cd $RR
   mkdir $GO_SAMPLE
   cd $GO_SAMPLE
   echo $(pwd)

   sbatch $SRCDIR_INI/02_STAR_run_query_STEP2.sh $AA/$RNAFILE_A \
                                                 $AA/$RNAFILE_B \
                                                 ${RNAFILE_A} \
                                                 ${RNAFILE_B} \
                                                 ${RR}/${GO_SAMPLE} \
                                                 ${GO_SAMPLE} \
                                                 ${REFG} \
                                                 $ST1



   sleep 1                                                  #pauses for 1 sec
   echo

done


squeue -u $USER                                    #check job status
echo

