#!/bin/bash

###### START UPPMAX RUNS

echo
echo "Starting Uppmax jobs ..."
echo





##### Paths and folders

# Path to folder containing input data
AA=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1e_cutadapt

# Path to output folder
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1f_condetri

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



### Number of libraries
N_LIBS=${#SAMPLE_LIST[@]}






######## file names


declare -A READS

for i in `seq 1 1 $N_LIBS`; do
   READS[$i,1]=${SAMPLE_LIST[$i]}_R1_001_cutadapt.fq.gz
   READS[$i,2]=${SAMPLE_LIST[$i]}_R2_001_cutadapt.fq.gz
done












# prepare output folder

#SRCDIR=$(pwd)                                    #remember path

#cd $PATH_MAIN/                                   #cd to project folder
#OUT_FOLDER=1e_QC_condetri2/                       #main output folder

#if [ -d "$OUT_FOLDER" ]                          #check whether main output folder already exists
#   then
#      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
#   else
#      mkdir $OUT_FOLDER                          #create output folder
#fi



# start jobs

cd $SRCDIR_INI

for i in `seq 1 1 $N_LIBS`; do                #loop starting all jobs
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   GO_SAMPLE=${SAMPLE_LIST[$i]}

   echo $RNAFILE_A
   echo $RNAFILE_B
   echo $GO_SAMPLE

   sbatch 01_condetri_run_query.sh $AA/$RNAFILE_A \
                                   $AA/$RNAFILE_B \
                                   ${RNAFILE_A} \
                                   ${RNAFILE_B} \
                                   $RR \
                                   ${GO_SAMPLE}

   sleep 1                                       #pauses for 1 sec
   echo
done

echo
squeue -u $USER                                    #check job status
echo
(END)

