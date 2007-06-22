#!usr/bin/perl -w

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

use Bio::SeqIO;                                            # Use Bio::SeqIO for the reading and writing of seq files in different formats
use Cwd;                                                   # Use the Cwd module to get the current working directory

# Get command-line arguments, or die with a usage statement
my $usage = "Fasta2Dir.pl infile infileformat outfile \n";
my $infile = shift or die $usage;                          # The full path of the infile that will be transformed
my $infileformat = shift or die $usage;                    # The format that the input file is in. Any valid bioperl format (ie. abi, ace, gcg, genbank, fasta, swiss, tigr etc.) 
my $outfile = shift or die $usage;                         # The name that will be used as the prefix for the output files. (ie. Output0000001, Output0000002, ..etc)
my $CurrentDir = cwd();                                    # Get the current working directory
my $OutputDir = $CurrentDir."/".$outfile."/";              # Set the output dirctory that the files will be created in
mkdir $OutputDir;                                          # Make the output directory, if the output directory exists the existing directory will not be overwritten

$SeqNum = 1;

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);


 # write each entry in the input file to the new output file
while (my $inseq = $seq_in->next_seq)
{

      $num = sprintf("%7d", $SeqNum);                      # Pad number with zeros so that the total length of the string is 7 characaters (ie. 0000012)
      $num=~ tr/ /0/;


      # Could possibly retrieve fasta file information here
      print "Primary ID: ".$inseq->primary_id."\n";
      my $info = $inseq->primary_id;
      @info_split = split(/\|/, $info);
      print "Parsed:".@info_split[3]."\n";

      my $seq_out = Bio::SeqIO->new('-file' => ">$OutputDir$outfile".@info_split[3],
                                    '-format' => fasta);
      print $outfile.$num." \n";                           # Print the number of the file so that we can keep track of where we are

      $seq_out->write_seq($inseq);                         # Write the individual fasta file

      $SeqNum++;                                           # Increment SeqNumber
}

print "The files have all been placed in the directory: \n";
print $OutputDir."\n";

exit;
