#!/bin/bash -l
#SBATCH -J Mapping_STAR_2
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 06:00:00
#SBATCH -A snic2019-3-35
#SBATCH -M snowy
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load star/2.5.3a
module load samtools/1.6




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

#folder containing results from STAR mapping STEP1
ST1=${8:?msg}






# Run STAR in two-pass mode







# Get list of SJ.out.tab files obtained for each sample during STEP 1

cd ${ST1}
find * -path "${SNAME:0:3}*" -maxdepth 1 -type f > $RR/list_all_files_OUT.txt       #gets names of all files present in output data folder(s) from STEP1, saves names to file

sjf="$(grep "\SJ.out.tab" $RR/list_all_files_OUT.txt)"                      #finds and concatenate paths associated to all SJ.out.tab files

echo
echo "SJ.out.tab files"
echo $sjf
echo





# SECOND PASS (to be performed, for each sample, after the first pass has been completed for all samples)
# (see section 8.1 in manual)


cd ${ST1}

STAR \
      --runMode alignReads \
      --runThreadN 10 \
      --genomeDir $REFG \
      --readFilesIn $READ1 $READ2 \
      --readFilesCommand zcat \
      --outFileNamePrefix $SNIC_TMP/ \
      --sjdbFileChrStartEnd $sjf \
      --outSAMattrIHstart 0 \
      --outFilterMultimapNmax 1 \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0

#      --outSAMstrandField intronMotif \













# Notes (STAR, mapping step)
#--runMode alignReads                           map reads
#--runThreadN N                                         specifies the number of threads that will be used by the program
#--outSAMstrandField intronMotif        adds information (to the SAM output file) required for downstream analysis with Cufflinks
#--genomeDir /path/to/index             specifies the directory containing the pre-built genome index
#--readFilesIn /path/to/reads/sample_1.fastq /path/to/reads/sample_2.fastq              is where you should list the FASTQ files that you wish to map
#--readFilesCommand zcat            uncompress gz files
#--outFileNamePrefix outDir             specifies the output directory
# --outSAMattrIHstart           default: 1      int>=0: start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie
#                                               The number of loci Nmap a read maps to is given by NH:i:Nmap field. Value of 1 corresponds to
#                                               unique mappers, while values >1 corresponds to multi-mappers. HI attrbiutes enumerates multiple
#                                               alignments of a read starting with 1 (this can be changed with the --outSAMattrIHstart - setting
#                                               it to 0 may be required for compatibility with downstream software such as Cufflinks or StringTie).
# --outFilterMultimapNmax               default: 10             int: maximum number of loci the read is allowed to map to. Alignments (all of
#                                                                                               them) will be output only if the read maps to no more loci than this value.
#                                                                                               Otherwise no alignments will be output, and the read will be counted as
#                                                                                               "mapped to too many loci" in the Log.final.out .


## >> The following should be set to zero if one gets a high percentage of reads being filtered out during mapping because such reads were deemed to be too short
##    (see also https://github.com/alexdobin/STAR/issues/169)
#--outFilterScoreMin                    default: 0              int: alignment will be output only if its score is higher than or equal to this value.
#--outFilterScoreMinOverLread   default: 0.66   real: same as outFilterScoreMin, but normalized to read length (sum of mates' lengths for paired-end reads)
#--outFilterMatchNmin                   default: 0              int: alignment will be output only if the number of matched bases is higher     than or equal to this value.
#--outFilterMatchNminOverLread  default: 0.66   real: sam as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads).








# Use Samtools to convert sam files into sorted bam format, and then to sort (by coordinate) and index the bam file
cd $SNIC_TMP
samtools view -bS Aligned.out.sam > Aligned.out.bam
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam
samtools index Aligned.out.sorted.bam


# Use Samtools to glean statistics about sorted BAM file
samtools flagstat Aligned.out.sorted.bam > STAR_flagstat.txt
#samtools idxstats Aligned.out.sorted.bam > STAR_idxstats.txt










##################### FILTERING >> NOT BEING DONE (read notes on #method 2 (STAR), below)




############ Filter out non-properly paired reads, reads that map more than once, and chimeric alignments
#                       Properly paired reads are assigned the 0x2 code by BWA
#           Bwa adds an XA auxiliary tag to reads that have another "valid" mapping (secondary)  and  "grep -v 'XA:Z:'" filters out all such reads.
#           "grep -v 'SA:Z:'" filters out supplementary alignments.
#           note: using '-F 0x100 -F 0x800' gets rid of secondary and supllementary reads, but not of primary reads that have secondary/supplementary reads.

#                       Definitions
#                       Properly paired: both reads of a pair map to the same chromosome/scaffold, oriented towards each other, and with a sensible insert size.
#                       Secondary: number of reads mapping to more than one position
#                       Supplementary: number of chimeric alignments (different parts of the read align to different parts of the genome)
#                       Duplicates: number of reads with identical copies present in library (duplicate reads will be removed later on)


# method 1 (when using BWA/BOWTIE2 aligners)
## if including only properly paired reads
#samtools view -hf 0x2 Aligned.out.sam \
#                  | grep -v -e 'XA:Z:' -e 'SA:Z:' \
#               > Aligned.out.filtered.sam

## if keeping non-properly paired reads
#samtools view -h $SNIC_TMP/${SNAME}.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $SNIC_TMP/${SNAME}.filtered.alignments.sam





# method 2 (STAR)  >>> not required as we filtered out multiloci reads while running STAR)
#cd $SNIC_TMP

# STAR flag for uniquely mapped reads:  HI:i:0
# STAR comment lines start with @
# NOTE:  Non-properly-paired reads are NOT being filtered out because in some cases we are using a reference genome from a sister species

#samtools view -h -F 0x100 -F 0x800 Aligned.out.sam \
#                  | grep -e 'HI:i:0' -e '^@' \
#               > Aligned.out.filtered.sam


#samtools view -h -F 0x100 -F 0x800 Aligned.out.sam \
#               > Aligned.out.filtered.sam





# Use Samtools to convert sam files into sorted bam format, and then to sort and index the bam file
#cd $SNIC_TMP
#samtools view -bS Aligned.out.filtered.sam > Aligned.out.filtered.bam
#samtools sort Aligned.out.filtered.bam -o Aligned.out.filtered.sorted.bam
#samtools index Aligned.out.filtered.sorted.bam


# Use Samtools to glean statistics about sorted BAM file
#samtools flagstat Aligned.out.filtered.sorted.bam > STAR_flagstat.filtered.txt
#samtools idxstats Aligned.out.sorted.bam > STAR_idxstats.txt








### copy files to results folder

# note: keep SAM file as it will be scrubbed downstream (removal of problematic reads)

cd $RR
rsync -ah $SNIC_TMP/Aligned.out.sorted.bam .
rsync -ah $SNIC_TMP/Aligned.out.sorted.bam.bai .
gzip -c $SNIC_TMP/Aligned.out.sam > $RR/Aligned.out.sam.gz
#rsync -ah $SNIC_TMP/Aligned.out.sam .
#rsync -ah $SNIC_TMP/Aligned.out.filtered.sorted.bam .
#rsync -ah $SNIC_TMP/Aligned.out.filtered.sorted.bam.bai .
rsync -ah $SNIC_TMP/Log* .
rsync -ah $SNIC_TMP/*flagstat* .
rsync -ah $SNIC_TMP/SJ.out.tab .




# runtime:  2h30m        (10 cores);
#           memory requirements > needs at least 6 cores
#           12h       (10 cores) if reference genome is Lepidium_meyenii (genome index was constructed differently as GFF file did not contain axons)
#      peak memory: 16GB > needs at least 3 cores


