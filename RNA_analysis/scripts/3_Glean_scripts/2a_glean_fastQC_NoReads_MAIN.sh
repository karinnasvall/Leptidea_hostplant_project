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
AA=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/1_fastqc

# Sample list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt

#results folder	
RR=$SRCDIR	

#output file
OUT_file=Total_Sequences_OUT1_Raw.txt									







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
   ZFILES[$i,1]=${SAMPLE_LIST[$i]}_L002_R1_001_fastqc.zip
   ZFILES[$i,2]=${SAMPLE_LIST[$i]}_L002_R2_001_fastqc.zip
done





##### Get data for each library pair


for i in `seq 1 1 $N_LIBS`; do                  #loop starting all jobs

   # unzip file
   zfile1=${ZFILES[$i,1]}
   zfile2=${ZFILES[$i,2]}
   echo
   echo $zfile1
   echo $zfile2
   cd $SRCDIR
   cp $AA/$zfile1 ./
   cp $AA/$zfile2 ./
   unzip -q $SRCDIR/$zfile1
   unzip -q $SRCDIR/$zfile2

   # cd to unzipped folder 1
   cd $SRCDIR/${zfile1%.zip}

   # glean filename
   filename_AUX1="$(grep "Filename" fastqc_data.txt | cut -f 2)"

   # glean total number of sequences
   tns_AUX1="$(grep "Total Sequences" fastqc_data.txt | cut -f 2)"

   # cd to unzipped folder 2
   cd $SRCDIR/${zfile2%.zip}

   # glean filename
   filename_AUX2="$(grep "Filename" fastqc_data.txt | cut -f 2)"

   # glean total number of sequences
   tns_AUX2="$(grep "Total Sequences" fastqc_data.txt | cut -f 2)"

   # store results
   RES_R1[i]="$filename_AUX1"
   RES_tns1[i]="$tns_AUX1"
   RES_tns2[i]="$tns_AUX2"

   # remove unzipped folder
   cd ..   
   rm -r -f $SRCDIR/${zfile1%.zip}
   rm -f $SRCDIR/$zfile1
   rm -r -f $SRCDIR/${zfile2%.zip}
   rm -f $SRCDIR/$zfile2
 
done





##### save to file


   
# put results in the right order and save to file
for i in `seq 1 1 $N_LIBS`; do
   if [[ "${i}" = "1" ]]
   then                  
       echo ${SAMPLE_LIST[$i]}_R1 ${RES_tns1[$i]} > $RR/$OUT_file
       echo ${SAMPLE_LIST[$i]}_R2 ${RES_tns2[$i]} >> $RR/$OUT_file
   else
       echo ${SAMPLE_LIST[$i]}_R1 ${RES_tns1[$i]} >> $RR/$OUT_file
       echo ${SAMPLE_LIST[$i]}_R2 ${RES_tns2[$i]} >> $RR/$OUT_file
   fi
done




#runtime:



