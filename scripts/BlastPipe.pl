#!/usr/bin/perl -w
# BLAST PIPELINE
# 04/11/2007
# J. Estill
# To run the Blast_PAWS program

# Directory of files to analyze
#my $AnalDir = "/home/jestill/projects/wheat_annotation/200704/tabblast/";
my $AnalDir = "/home/jestill/projects/wheat_annotation/20070619/bac2blast/";

#system("./Blast_PAWS.pl -i $AnalDir -d embl -g");
#system("./Blast_PAWS.pl -i $AnalDir -d estn -g");
#system("./Blast_PAWS.pl -i $AnalDir -d estx -g");
#system("./Blast_PAWS.pl -i $AnalDir -d prot -g");
#system("./Blast_PAWS.pl -i $AnalDir -d repn -g");
#system("./Blast_PAWS.pl -i $AnalDir -d repx -g");
#system("./Blast_PAWS.pl -i $AnalDir -d cont -g");
#system("./Blast_PAWS.pl -i $AnalDir -d genx -g");


# 06/25/2007
# The following is a temporary fix to do proper parsing 
# of the GFF files for BLAST results
system("./GFF_Blast_PAWS.pl -i $AnalDir -d embl -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d estn -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d estx -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d prot -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d repn -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d repx -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d cont -g");
system("./GFF_Blast_PAWS.pl -i $AnalDir -d genx -g");

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 04/11/2007
# - Set up pipeline to analyze a set of wheat BACs
