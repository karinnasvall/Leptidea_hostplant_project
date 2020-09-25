#!/usr/bin/perl -w

#use warnings;
#use strict;

# Add gene names from STAR+stringtie (filtered) outfile to outfile from DEseq (filtered for SBG) and print to new file
# Used for making volcano plot with gene names for sex biased genes

# Usage: perl 01_add_gene_names.pl DEseq.file STAR.file out.file

# Lars Höök 2017


$DEseq=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2/04_SpaSwe_Hostplant_effect/Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_SpaSwe_Host.txt/
$STAR=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/5_Get_Gene_Names_scripts/STAR_Aligned_out_Gene_Names_list.tab/
$outfile=/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/5_Get_Gene_Names/01_SpaSwe_Host/



if ( @ARGV != 3 ) {
    print "Usage: perl add_gene_names_to_SBG.pl DEseq.file STAR.file out.file\n";
    exit;
}


my ($DEseq, $STAR, $outfile) = @ARGV ;

open (DES, $DEseq) or die "Can't open $DEseq\n";
while (my $line = <DES>) {

        chomp $line;
        $line =~ s/\Q"//g;      #remove ""
        my @col_1 = split(/\,/, $line); #split by comma
        my $gene_id1 = $col_1[0];
        my $baseMean = $col_1[1];

                if ($baseMean eq "baseMean") {

                open(OUT, ">>$outfile");
                print OUT "$line", ",", "GeneName\n";

                }

                        else {

                        open (STA, $STAR) or die "Can't open $STAR\n";
                        while (my $genes = <STA>) {             #read through DEseq file

                        my @col_2 = split(/\t/, $genes); #split by tab and pick gene id, lfc and p-adj
                        my $gene_id2 = $col_2[0];
                        my $gene_name = $col_2[1];
                        $gene_name =~ s/\Q-/ /g;
                        $gene_name =~ s/\Q./ /g;

                                if ($gene_id1 eq $gene_id2) {   #print line to outfile if gene_id matches

                                my $line_name = join ",", $line, $gene_name;
                                open(OUT, ">>$outfile");
                                print OUT "$line_name\n";

                                }

                                else { next }

                        }

                }

}

