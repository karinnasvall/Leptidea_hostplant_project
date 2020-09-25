#!/bin/bash -l
#SBATCH -J STAR-index
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 10:00:00
#SBATCH -A snic2019-3-35
##SBATCH -M snowy
#SBATCH --mail-user veronika.mrazek@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load star/2.5.3a






##### Paths and folders

# Path to folder containing reference genome to be indexed
AA=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/LS_genome

# Reference genome to be indexed
RG=NB_Leptidea.filtered.fa

# Path to annotation file
ANN=/crex/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/LS_genome/leptidea_sinapis_rc1.gff

#remember initial path
SRCDIR_INI=$(pwd)





echo
echo "Genome assemply:" ${AA}/${RG}
echo "Annotation file:" ${ANN}
echo


cd ${AA}

mkdir star_index

STAR \
    --runMode genomeGenerate \
    --runThreadN 3 \
    --genomeDir ${AA}/star_index/ \
    --genomeFastaFiles ${AA}/${RG} \
    --sjdbGTFfile ${ANN} \
    --sjdbGTFtagExonParentTranscript Parent



    # --sjdbGTFfeatureExon CDS              >> use when indexing Lepidium_Meyenii (see note A)

    # --genomeChrBinNbits 7                     >> use when indexing Lepidium_Meyenii (see note B)

    ## Note A: Lepidium_Meyenii annotation file contains no exon marks (unlike all the other gff3 files associated to the other genomes). Instead, we will use CDSs for building transcripts.
    ## Note B: If you are using a genome with a large > 5,000 number of references (chrosomes/scaffolds), you may need to reduce the --{genomeChrBinNbits to reduce RAM consumption.
    ## The following scaling is recommended: --genomeChrBinNbits} = min(18, log2(GenomeLength/NumberOfReferences)).
    ## For example, for 3~gigaBase genome with 100,000 chromosomes/scaffolds, this is equal to 15.









#Notes (STAR genomeGenerate)
#--runMode                      genomeGenerate option directs STAR to run genome indices generation job.
#--runThreadN N         specifies the number of threads that will be used by the program
#--genomeDir            specifies path to the directory (henceforth called "genome directory" where the
#                                       genome indices are stored. This directory has to be created (with mkdir) before STAR run
#                                       and needs to writing permissions.
#--genomeFastaFiles specified one or more FASTA files with the genome reference sequences
#--sjdbOverhang         specifies the length of the genomic sequence around the annotated junction
#                                       to be used in constructing the splice junctions database. Ideally, this length should be equal
#                                       to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
#                                       2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
#                                       ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
#                                       well as the ideal value. [USE ONLY when annotation file is also provide]
# --sjdbGTFfile         specifies the path to the file with annotated transcripts in the standard GTF
#                                       format. STAR will extract splice junctions from this file and use them to greatly improve
#                                       accuracy of the mapping. While this is optional, and STAR can be run without annotations,
#                                       using annotations is highly recommended whenever they are available.
# --sjdbOverhang        specifies the length of the genomic sequence around the annotated junction
#                                       to be used in constructing the splice junctions database. Ideally, this length should be equal
#                                       to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
#                                       2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
#                                       ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
#                                       well as the ideal value.

# Annotations in GFF format.
# In addition to the aforementioned options, for GFF3 formatted annotations you need to use
# --sjdbGTFtagExonParentTranscript Parent. In general, for --sjdbGTFfile files STAR only
# processes lines which have --sjdbGTFfeatureExon (=exon by default) in the 3rd field (column). The exons are assigned to the
# transcripts using parent-child relationship defined by the
# --sjdbGTFtagExonParentTranscript (=transcript id by default) GTF/GFF attribute

