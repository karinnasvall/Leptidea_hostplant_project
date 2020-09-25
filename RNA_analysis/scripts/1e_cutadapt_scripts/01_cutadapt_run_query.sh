#!/bin/bash -l
#SBATCH -J Cutadapt
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -A snic2019-3-35
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL




# load modules
module load bioinfo-tools
module load cutadapt/1.16
module load FastQC/0.11.5
#module load MultiQC/1.6.   >>> not compatible with cutadapt (uses python 2, while cutadapt uses python 3)





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
#SNAME=${SNAME}_C1




##unzip fq.gz files
gunzip -c $READ1>$SNIC_TMP/LIB1.fastq
gunzip -c $READ2>$SNIC_TMP/LIB2.fastq

#cd $SNIC_TMP                                # use local scratch disk





##### run cutadapt


# remove Poly-A missed by prinseq (such as internal Poly-A/N or Poly-A/N with errors)

cutadapt \
               --quality-base=33 \
           --quality-cutoff=30,30 \
           -b "A{10}" \
           -b "T{10}" \
           -B "T{10}" \
           -B "A{10}" \
           --times=2 \
           --overlap=3 \
           --trim-n \
           --minimum-length=30 \
           --max-n=0.1 \
           --nextseq-trim=30 \
               --output $SNIC_TMP/${SNAME}_R1_001_cutadapt.fq.gz \
                   --paired-output $SNIC_TMP/${SNAME}_R2_001_cutadapt.fq.gz \
                   $READ1 \
                   $READ2


# trim 2 nucleotide at 5' end (FASTQC analysis shows unbalanced nucleotide composition at the end of read)

#cutadapt \
#               --quality-base=33 \
#               -u -2 \
#                   -U -2 \
#           --quality-cutoff=30,30 \
#           --trim-n \
#           --minimum-length=30 \
#           --max-n=0.1 \
#           --pair-filter=any \
#               --output $SNIC_TMP/${SNAME}_R1_001_cutadapt.fq.gz \
#                   --paired-output $SNIC_TMP/${SNAME}_R2_001_cutadapt.fq.gz \
#                   $SNIC_TMP/${SNAME}_R1_001_cutadapt_AUX.fq.gz \
#                  $SNIC_TMP/${SNAME}_R2_001_cutadapt_AUX.fq.gz



#           -b "A{100}" \
#           -b "T{100}" \
#           -B "T{100}" \
#           -B "A{100}" \




### copy results to output folder
cd $RR
cp $SNIC_TMP/${SNAME}_R1_001_cutadapt.fq.gz .
cp $SNIC_TMP/${SNAME}_R2_001_cutadapt.fq.gz .




### perform FastQC analysis
fastqc -o $RR/ ${SNAME}_R*_001_cutadapt.fq.gz


# create single report
#multiqc $RR/ --cl_config "max_table_rows: 10000"











### Notes:
#     -u LENGTH, --cut=LENGTH                           Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the
#                                                               beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different
#                                                               signs. This is applied *before* adapter trimming.
#     -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff=[5'CUTOFF,]3'CUTOFF                      Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if
#                                                                                                                               data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given,
#                                                                                                                               the 5' end is trimmed with the first cutoff, the 3' end with the second.
#     --quality-base=QUALITY_BASE                       Assume that quality values in FASTQ are encoded as ascii(quality + QUALITY_BASE).
#                                                                                       This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
#       -a ADAPTER, --adapter=ADAPTER                   Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent
#                                                               bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.
#   -g ADAPTER, --front=ADAPTER                         Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any
#                                                               preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended
#     -O MINLENGTH, --overlap=MINLENGTH         Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
#                                                               ('anchoring'), the adapter is only found if it is a prefix of the read.
#     -n COUNT, --times=COUNT                           Remove up to COUNT adapters from each read. Default: 1
#     --trim-n                              Trim N's on ends of reads.
#     --max-n=COUNT                                             Discard reads with more than COUNT 'N' bases.
#                                                                                       If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
#     -m LENGTH, --minimum-length=LENGTH        Discard reads shorter than LENGTH. Default: 0
#     --minimum-length LENGTH1:LENGTH2 or -m LENGTH1:LENGTH2                    When trimming paired-end reads, the minimum lengths for R1 and R2 can be specified separately by separating them with a colon (:).
#                                                                                                                                               If the colon syntax is not used, the same minimum length applies to both reads, as discussed above. Also, one of the values can be
#                                                                                                                                               omitted to impose no restrictions. For example, with -m 17:, the length of R1 must be at least 17, but the length of R2 is ignored.
#     --pair-filter=(any|both)                          Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
#     -o FILE, --output=FILE                            Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent
#                                                       to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
#     -p FILE, --paired-output=FILE                     Write second read in a pair to FILE.



# Each read is processed in the following order:

# 1.    Read modification options are applied. This includes adapter removal,  quality trimming,  read name modifications etc.
#               The order in which they are applied is the order in which they are listed in the help shown by cutadapt --help under the “Additional read modifications” heading.
#               Adapter trimming itself does not appear in that list and is done after quality trimming and before length trimming (--length/-l).
#
#               Additional read modifications:
#                               -u LENGTH, --cut=LENGTH
#                               --nextseq-trim=3'CUTOFF
#                               -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff=[5'CUTOFF,]3'CUTOFF
#                               --quality-base=QUALITY_BASE
#                               -l LENGTH, --length=LENGTH
#                               --trim-n
#                               --length-tag=TAG
#                               --strip-suffix=STRIP_SUFFIX
#                               -x PREFIX, --prefix=PREFIX

# 2.    Filtering options are applied, such as removal of too short or untrimmed reads. Some of the filters also allow to redirect a read to a separate output file.
#               The filters are applied in the order in which they are listed in the help shown by cutadapt --help under the “Filtering of processed reads” heading.
#
#               Filtering of processed reads:
#                               -m LENGTH, --minimum-length=LENGTH
#                               -M LENGTH, --maximum-length=LENGTH
#                               --max-n=COUNT
#                               --discard-trimmed, --discard
#                               --discard-untrimmed, --trimmed-only
#
# 3.    If the read has passed all the filters, it is written to the output file.












#### Runtime: 3h




