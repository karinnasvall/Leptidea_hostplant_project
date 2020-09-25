#!/bin/bash -l
#SBATCH -J fastQC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2019-3-35
##SBATCH -M snowy
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
##SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=ALL



# load modules
module load bioinfo-tools
module load FastQC/0.11.5
module load MultiQC/1.6



## PATHS and filenames

#initial folder
SRCDIR=$(pwd)                                                   

#read name of files containing paired reads (with path)
READ1=${1:?msg}		
READ2=${2:?msg}

#read name of files containing paired reads (without path)
rawReads1=${3:?msg}     
rawReads2=${4:?msg}

echo $rawReads1
echo $rawReads2


# output folder (root name)
RR=${5:?msg}

#sample name
SNAME=${6:?msg}




## Run FastQC

fastqc -o $RR/ $READ1
fastqc -o $RR/ $READ2





# create single report
#multiqc ./ --cl_config "max_table_rows: 10000"


# runtime (for 2 fasta.gz file): 16 min



