#!/bin/bash -l
#SBATCH -J fastq_masker
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2019-3-35
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load Fastx/0.0.14
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







##### unzip fq.gz files to scratch disk
gunzip -c $READ1>$SNIC_TMP/READ1.fastq
gunzip -c $READ2>$SNIC_TMP/READ2.fastq





##### use fastq_masker to mask low quality nucleotides (mask nucleotides with quality below 20)
##### each file (paired end) masked independently

fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $SNIC_TMP/READ1.fastq \
        -o $SNIC_TMP/${SNAME}_R1_001_masked.fq.gz


fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $SNIC_TMP/READ2.fastq \
        -o $SNIC_TMP/${SNAME}_R2_001_masked.fq.gz






# notes:
#   [-q N]       = Quality threshold - nucleotides with lower quality will be masked
#                  Default is 10.
#   [-r C]       = Replace low-quality nucleotides with character C. Default is 'N'
#   [-z]         = Compress output with GZIP.
#   [-i INFILE]  = FASTQ input file. default is STDIN.
#   [-o OUTFILE] = FASTQ output file. default is STDOUT.
#   [-v]         = Verbose - report number of sequences.
#                  If [-o] is specified,  report will be printed to STDOUT.
#                  If [-o] is not specified (and output goes to STDOUT),
#                  report will be printed to STDERR.







### copy masked libraries to results folder

cd $RR
cp $SNIC_TMP/${SNAME}_R1_001_masked* .
cp $SNIC_TMP/${SNAME}_R2_001_masked* .




### perform FastQC analysis
fastqc -o $RR/ $RR/${SNAME}*




# create single report
#multiqc $RR/ --cl_config "max_table_rows: 10000"


# runtime: 2h

