#!/usr/bin/perl -w
#-----------------------------------------------------------+
# Fasta Extractor
#
#-----------------------------------------------------------+

# Author: James C. Estill
# Contact: jestill@uga.edu
# Started: 7/22/04
# Last Visited: 02/20/2007
#
# Takes an input file in any valid bioperl format and converts it
# to a directory containing one fasta file for each sequence. This was
# written to create a set of files for use on the cluster. Files are currently
# placed in the directory that the script is run under.
#

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;
use Cwd;

#-----------------------------+
# GET VARS FROM COMMAND LINE  |
#-----------------------------+
my $usage = "x2fasta.pl infile infileformat outfile \n";
my $infile = shift or die $usage;
my $infileformat = shift or die $usage;
my $outfile = shift or die $usage;

my $CurrentDir = cwd();
my $infile = $CurrentDir."/".$infile;
my $outfile = $CurrentDir."/".$outfile;

print $CurrentDir."\n";
print "IN".$infile."\n";
print "OUT".$outfile."\n";
#exit;

my $SeqNum = 0;

#-----------------------------+
# OPEN I/O FILES              |
#-----------------------------+
# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat) ||
    die "Can not open input file:n\$infile\n";

my $seq_out = Bio::SeqIO->new('-file' => ">$outfile",
			      '-format' => 'fasta') ||
    die "Can not open output file:\n$outfile\n";

#-----------------------------+
# CONVERT INPUT FILE FORMAT   |
# TO OUTPUT FILE FORMAT       |
#-----------------------------+
while (my $inseq = $seq_in->next_seq)
{
      $SeqNum++;
      print "Converting seq: ".$SeqNum."\n";
      $seq_out->write_seq($inseq);
}

exit;

#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 02/20/2007
# - Updated program to my current coding conventions
#   and basically rewrote the entire program
