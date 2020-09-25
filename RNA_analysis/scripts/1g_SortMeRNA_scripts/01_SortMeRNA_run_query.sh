#!/bin/bash -l
#SBATCH -J SortMeRNA
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 04:00:00
#SBATCH -A snic2019-3-35
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL




# load modules
module load bioinfo-tools
module load SortMeRNA/2.1b
#module load bowtie2/2.3.4.1



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

## path to SortMeRNA scripts
SMRS=/sw/apps/bioinfo/SortMeRNA/2.1b/rackham/scripts/





# unzip files
gunzip -c $READ1 > $SNIC_TMP/${rawReads1%.gz}
gunzip -c $READ2 > $SNIC_TMP/${rawReads2%.gz}



# merge forward and reverse files
$SMRS/merge-paired-reads.sh $SNIC_TMP/${rawReads1%.gz} $SNIC_TMP/${rawReads2%.gz} $SNIC_TMP/SMRNA_IN.fastq





##### filter reads
##### note: the 'Total reads passing E-value threshol' value in the log output file will count double for the same read if it maps to two screening labraries (count three times if it maps to 3 screening libraries, etc)

sortmerna -a 10\
          -m 50000\
          --ref $SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90:\
$SORTMERNA_DBS/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DBS/index/silva-bac-23s-id98:\
$SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95:\
$SORTMERNA_DBS/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DBS/index/silva-arc-23s-id98:\
$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95:\
$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98:\
$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98:\
$SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98\
          --reads $SNIC_TMP/SMRNA_IN.fastq\
          --paired_in\
          --blast '1'\
          --num_alignments 1\
          --fastx\
          --aligned $SNIC_TMP/${SNAME}_SORTMERNA_rRNA\
          --other $SNIC_TMP/SMRNA_OUT\
          --log\
          -v





## unmerge paired reads back into two separate files
$SMRS/unmerge-paired-reads.sh $SNIC_TMP/SMRNA_OUT.fastq $SNIC_TMP/${SNAME}_R1_001_SORTMERNA_CLEAN.fastq $SNIC_TMP/${SNAME}_R2_001_SORTMERNA_CLEAN.fastq







## copy filtered libraries to results folder
cd $RR
cp $SNIC_TMP/*.log .
cp $SNIC_TMP/*.blast .
gzip -c $SNIC_TMP/${SNAME}_R1_001_SORTMERNA_CLEAN.fastq > $RR/${SNAME}_R1_001_SORTMERNA_CLEAN.fastq.gz
gzip -c $SNIC_TMP/${SNAME}_R2_001_SORTMERNA_CLEAN.fastq > $RR/${SNAME}_R2_001_SORTMERNA_CLEAN.fastq.gz
cp $SNIC_TMP/*SORTMERNA_rRNA* .






# notes (sortmerna):
# -a                INT             number of threads to use
# -m                INT             INT Mbytes for loading the reads into memory (default: 1024 , maximum: 128693)
# --ref             STRING,STRING   FASTA reference file, index file
#                                         (ex. --ref /path/to/file1.fasta,/path/to/index1)
#                                         If passing multiple reference files, separate
#                                         them using the delimiter ':',
#                                         (ex. --ref /path/to/file1.fasta,/path/to/index1:/path/to/file2.fasta,path/to/index2)
# --fastx           BOOL            output FASTA/FASTQ file (for aligned and/or rejected reads)
# --num_alignments  INT             report first INT alignments per read reaching E-value (--num_alignments 0 signifies all alignments will be output)
# --best            INT             report INT best alignments per read reaching E-value by searching --min_lis INT candidate alignments
#                                       (--best 0 signifies all candidate alignments will be searched)
# --blast           STRING          output alignments in various Blast-like formats
#                                        '0' - pairwise
#                                        '1' - tabular (Blast -m 8 format)
#                                        '1 cigar' - tabular + column for CIGAR
#                                        '1 cigar qcov' - tabular + columns for CIGAR
#                                                         and query coverage
#                                        '1 cigar qcov qstrand' - tabular + columns for CIGAR,
#                                                                query coverage and strand
#
# --otu_map         BOOL            output OTU map (input to QIIME's make_otu_table.py)            off
# --log             BOOL            output overall statistics
# -v                BOOL            verbose
# --aligned         STRING          aligned reads filepath + base file name (appropriate extension will be added)
# --other           STRING          rejected reads filepath + base file name (appropriate extension will be added)
# For users who wish to keep the order of their paired-ended reads, two options are available.  If one read aligns and the other one not then,
#                                       (1) --paired-in         will put both reads into the file specified by --aligned
#                                       (2) --paired-out        will put both reads into the file specified by --other
#                               The  first  option, --paired-in, is  optimal  for  users  that  want  all  reads  in  the --other file  to  be non-rRNA.
#                               However, there are small chances that reads which are non-rRNA will also be put into the --aligned file.
# --paired_in       BOOL            both paired-end reads go in --aligned fasta/q file             off
#                                         (interleaved reads only, see Section 4.2.4 of User Manual)
# --paired_out      BOOL            both paired-end reads go in --other fasta/q file               off
#                                         (interleaved reads only, see Section 4.2.4 of User Manual)
# --reads           STRING          FASTA/FASTQ reads file



# runtime: 2h

