#!/bin/bash -l
#SBATCH -J Mapping_STAR_1
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 06:00:00
#SBATCH -A snic2019-3-35
#SBATCH -M snowy
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load star/2.5.3a






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
echo $SRCDIR

# output folder (root name)
RR=${5:?msg}

#sample name
SNAME=${6:?msg}

#reference genome
REFG=${7:?msg}

echo
echo "reference genome"
echo $REFG

#Annotation file
#ANNF=${8:?msg}

#echo $ANNF










# Run STAR in two-pass mode

# For the most sensitive novel junction discovery, run STAR in the 2-pass mode. It does not increase the number
# of detected novel junctions, but allows to detect more splices reads mapping to novel junctions.
# The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions
# detected in the first pass, and use them as "annotated" junctions for the 2nd pass mapping.



# FIRST PASS (run for each sample individually)

cd $RR

STAR \
      --runMode alignReads \
      --runThreadN 3 \
      --outSAMstrandField intronMotif \
      --genomeDir $REFG \
      --readFilesIn $READ1 $READ2 \
      --readFilesCommand zcat \
      --outFileNamePrefix $RR/ \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0


rm -f $RR/Aligned.out.sam                       #file not required (final sam file produced during the second pass)

# Notes (STAR, mapping step)
#--runMode alignReads                           map reads
#--runThreadN N                                         specifies the number of threads that will be used by the program
#--outSAMstrandField intronMotif        adds information (to the SAM output file) required for downstream analysis with Cufflinks
#--genomeDir /path/to/index             specifies the directory containing the pre-built genome index
#--readFilesIn /path/to/reads/sample_1.fastq /path/to/reads/sample_2.fastq              is where you should list the FASTQ files that you wish to map
#--readFilesCommand zcat            uncompress gz files
#--outFileNamePrefix outDir             specifies the output directory

## >> The following should be set to zero if one gets a high percentage of reads being filtered out during mapping because such reads were deemed to be too short
##    (see also https://github.com/alexdobin/STAR/issues/169)
#--outFilterScoreMin                    default: 0              int: alignment will be output only if its score is higher than or equal to this value.
#--outFilterScoreMinOverLread   default: 0.66   real: same as outFilterScoreMin, but normalized to read length (sum of mates' lengths for paired-end reads)
#--outFilterMatchNmin                   default: 0              int: alignment will be output only if the number of matched bases is higher     than or equal to this value.
#--outFilterMatchNminOverLread  default: 0.66   real: sam as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads).





#note: runtime: < 3h  (3 cores)
#               12h (10 cores) if reference genome is Lepidium_meyenii (genome index was constructed differently as GFF file did not contain axons)
#      peak memory: 16GB > needs at least 3 cores


