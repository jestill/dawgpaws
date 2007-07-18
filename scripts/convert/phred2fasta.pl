#!/usr/bin/perl -w

# -----------------------------------------------------------------------
# Fasta2Dir.pl
#
# Author: James C. Estill
# Contact: jestill@uga.edu
# Started: 7/22/04
# Last Visited: 3/7/06
# Tried to get this to get the NCBI information from the fasta file.
#
# Takes an input file in any valid bioperl format and converts it
# to a directory containing one fasta file for each sequence. This was
# written to create a set of files for use on the cluster. Files are placed
# in a subdirectory under the current working directory that is given the
# same name as the outfile prefix.
#
# Example Usage:
#
# Fasta2Dir.pl /home/jestill/BigFastaFile.fasta fasta IndFiles.fasta
# Fasta2Dir.pl /home/jestill/ArabidopsisGenes.xml tigr ArabGenes.fasta
#
#------------------------------------------------------------------------

print "Phred to Fasta Coversion program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;
use Getopt::Std;
use Cwd;

#-----------------------------+
# Get command-line arguments, |
#-----------------------------+
my %Options;
getopts('i:o:', \%Option);      # Load the options that are sent into the hash named options

# GET START TIME

my $usage = "Usage is: phred2fasta.pl -i indir -o outfile \n";
my $indir = $Option{i} ||
    die "ERROR!!\nInput directory not specified\n$usage";  
my $outfile = $Option{o} ||
    die "ERROR!!\nOutput file not specified\n$usage";
my $CurrentDir = cwd();

#-----------------------------+
# GET LIST OF FILES TO COVERT |
#-----------------------------+
opendir( DIR, $indir );
my @fileList = grep !/^\.\.?$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# INITITIALIZE VARIABLES      |
#-----------------------------+
$SeqNum = 1;

foreach my $file ( @fileList )
{
    $SeqNum++;

    # Debug statement, just do first 10 records
    #if ($SeqNum == 10){exit;}
    print "Processing ".$SeqNum."\n";

    # OPEN INPUT FILE
    $FullPath = $indir.$file;
    my $seq_in = Bio::SeqIO->new('-file' => "<$FullPath",
				 '-format' => 'phd') ||
				     die "Can not open file".$FullPath;

    # COVERT THE INPUT SEQ TO FASTA AND APPEND TO OUTFILE
    while (my $inseq = $seq_in->next_seq)
    {

	my $seq_out = Bio::SeqIO->new('-file' => ">>$outfile",
				      '-format' => 'fasta');

	$seq_out->write_seq($inseq);                         # Write the individual fasta file

    }

}

print "The $SeqNum fileles have all been placed in the file:\n$outfile\n";

exit;
