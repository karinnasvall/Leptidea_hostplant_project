#!/bin/bash -l
#SBATCH -J condetri
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2019-3-35
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load FastQC/0.11.5






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




### Path to condetri
CONDETRI=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/1f_condetri_scripts/condetri_v2.3.pl


#unzip files
gunzip -c $READ1 > $SNIC_TMP/${rawReads1%.gz}
gunzip -c $READ2 > $SNIC_TMP/${rawReads2%.gz}










##### Run condetri
##### Usage:   perl condetri.pl -fastq1=file1 [-fastq2=file2 -prefix=s -cutfirst=i -cutlast=i -rmN -notrim -hq=i -lq=i -frac=[0,1] -lfrac=[0,1] -minlen=i -mh=i -ml=i -sc=i -pb=s]


cd $SNIC_TMP

PRXOUT=condetri-${SNAME}              #output file(s) prefix

perl $CONDETRI \
            -fastq1=$SNIC_TMP/${rawReads1%.gz} \
            -fastq2=$SNIC_TMP/${rawReads2%.gz} \
            -prefix=$PRXOUT \
            -cutfirst=0 \
            -cutlast=0 \
            -rmN \
            -hq=30 \
            -lq=0 \
            -frac=0.8 \
            -minlen=30 \
            -mh=1 \
            -ml=1 \
            -sc=33





### Compress filtered libraries and copy to output folder
gzip -c $SNIC_TMP/${PRXOUT}_trim1.fastq > $RR/${SNAME}_R1_001_condetri.fq.gz
gzip -c $SNIC_TMP/${PRXOUT}_trim2.fastq > $RR/${SNAME}_R2_001_condetri.fq.gz


### Run FASTQC
fastqc -o $RR/ $RR/${SNAME}_R*_001_condetri.fq.gz












# notes:
# -fastq1               first file of paired-end reads
# -fastq2               second file of paired-end reads
# -prefix               string Prefix for the output file(s)
# -cutfirst=i   Remove i bases from 5'-end before any trimming
# -cutlast=i    Remove i bases from the 3'-end before any trimming
# -rmN                  Remove Ns from 5' end before any trimming
# -hq:                  Bases are removed from the 3'-end if the quality score is lower
#                               than some hq quality threshold;
# -ml                   Max number of lq bases allowed within a stretch of hq bases
# -lq                   Low quality threshold (read filtered out if there is a base in the
#                               read with quality<lq)
# -mh                   When i consecutive hq bases is reached, the trimming stops
# -frac=[0,1]   Fraction of read that must exceed hq after quality trimming
# -lfrac=[0,1]  Maximum fraction of bases with qual<lq after quality trimming[0]
# -minlen       Min allowed read length
# -sc=i                 Illumina scoring table, Score=ASCII-sc. Used to be 64 for Illumina/Solexa.
#                               33 is Sanger standard, and also used for newer Illumina data (1.8+)





# create single report
#multiqc $RR/ --cl_config "max_table_rows: 10000"



# runtime: 5h



(END)

