#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# Fasta2Dir.pl                                              |
#  $Id:: Fasta2Dir.pl 63 2007-07-12 19:41:53Z JamesEstill $:|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/22/2004                                       |
# UPDATED: 07/11/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Takes an input file in any valid bioperl format and      |
#  converts it to a directory containing one fasta file for |
#  each sequence. Files are placed in a subdirectory under  |
#  the current working directory that is given the same     |
#  name as the outfile prefix.                              |
#                                                           |
# Usage:                                                    |
#                                                           |
#  Fasta2Dir.pl BigFastaFile.fasta fasta IndFiles.fasta     |
#  Fasta2Dir.pl ArabidopsisGenes.xml tigr InGenes.fasta     |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;               # Read and write seq files
use Cwd;                      # Get the current working directory

#-----------------------------+
# CMD LINE                    |
#-----------------------------+
my $usage = "Fasta2Dir.pl infile infileformat outfile \n";

my $infile = shift or die $usage;
my $infileformat = shift or die $usage;
my $outfile = shift or die $usage;  # The name that will be used 
                                    # as the prefix for the output files. 
                                    # (ie. Output0000001, Output0000002, ..etc)

# Get the current working directory
my $CurrentDir = cwd();
# Set the output dirctory that the files will be created in
my $OutputDir = $CurrentDir."/".$outfile."/";
mkdir $OutputDir;

$SeqNum = 1;

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);


 # write each entry in the input file to the new output file
while (my $inseq = $seq_in->next_seq) {
    
    # Pad number with zeros so that the total length of the 
    # string is 7 characaters (ie. 0000012)
    $num = sprintf("%7d", $SeqNum);
    $num=~ tr/ /0/;
    
    my $seq_out = Bio::SeqIO->new('-file' => ">$OutputDir$outfile".$num,
				  '-format' => fasta);
    
    print $outfile.$num." \n";
    
    $seq_out->write_seq($inseq);  # Write the individual fasta file  
    $SeqNum++;                    # Increment SeqNumber
}

print "The files have all been placed in the directory: \n";
print $OutputDir."\n";

exit;


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 07/11/2007
#  - Updated code for jperl
