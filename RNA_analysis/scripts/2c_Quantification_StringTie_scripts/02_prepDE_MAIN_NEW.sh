module load python/2.7.6

#OUTPUT FOLDER
#RR=/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/NEW_03_Normalize_counts_with_DESEQ2
RR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2c_Quantification_StringTie_STEP2

#remember current path
SRCDIR=$(pwd)



#### Run script
#    Note: for instructions on how to run script, see http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq

cd $RR

python $SRCDIR/prepDE.py \
                           --input=$SRCDIR/prepDE_samplelist.txt \
                           --length=262

# read length set to 2x131bp (where 131bp is the average read length for each library in paired-end mode)

# Notes:
# Options:
#-i INPUT, --input=INPUT, --in=INPUT    the parent directory of the sample sub-directories or a .txt file listing sample IDs and the paths to GTF files in
#                                       tab-delimited format [default: ballgown]
#-g G                                   where to output the gene count matrix [default: gene_count_matrix.csv]
#-t T                                   where to output the transcript count matrix [default: transcript_count_matrix.csv]
#-l LENGTH, --length=LENGTH             the average read length [default: 75]
#-p PATTERN, --pattern=PATTERN          a regular expression that selects the sample subdirectories
#-c, --cluster                          whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)
#-s STRING, --string=STRING             if a different prefix is used for geneIDs assigned by StringTie [default: MSTRG]
#-k KEY, --key=KEY                      if clustering, what prefix to use for geneIDs assigned by this script [default: prepG]
#--legend=LEGEND                        if clustering, where to output the legend file mapping transcripts to assigned geneIDs [default: legend.csv]



