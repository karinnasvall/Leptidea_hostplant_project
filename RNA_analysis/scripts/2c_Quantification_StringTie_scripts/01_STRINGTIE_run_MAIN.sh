###### START UPPMAX RUNS

# Performs transcript quantification with StringTie

echo
echo "Starting Uppmax jobs ..."
echo





##### Paths and folders

# Path to folder containing input data
AA=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2b_Mapping_STAR_STEP2

# Path to output folder
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2c_Quantification_StringTie_STEP1

# Sample list
SPL=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1a_QC_fastQC_scripts/sample_list_P11917.txt

# Annotation file
ANNF=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/LS_genome/leptidea_sinapis_rc1.gff



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














#start jobs

for i in `seq 1 1 $N_SAMPLES`; do

   BAMFILE_A=Aligned.out.sorted.bam
   #BAMFILE_A=Aligned.out.sorted.clean-2.bam
   GO_SAMPLE=${SAMPLE_LIST[$i]}
   GO_ANNOT=${ANNF}

   echo $AA/$GO_SAMPLE/$BAMFILE_A
   echo $GO_SAMPLE
   echo $GO_ANNOT

   cd $RR
   mkdir $GO_SAMPLE
   cd $GO_SAMPLE
   echo $(pwd)

   sbatch $SRCDIR_INI/01_STRINGTIE_run_query.sh $AA/$GO_SAMPLE/$BAMFILE_A \
                                                $BAMFILE_A \
                                                $RR/$GO_SAMPLE \
                                                $GO_SAMPLE \
                                                $GO_ANNOT


   sleep 1                                                  #pauses for 1 sec
   echo

done


squeue -u $USER                                    #check job status
echo

