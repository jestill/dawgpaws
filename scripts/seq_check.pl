#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq_check.pl - check sequence parameters                  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 03/09/2006                                       |
# UPDATED: 07/11/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Check some of the properties of the sequences files in   |
#  a valid biopel file format and report sequences that are |
#  outside of expected parameters.                          |
#                                                           |
# USAGE:                                                    |
# seq_check.pl /home/jestill/BigFastaFile.fasta fasta       |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;               # Read and write seq files in different formats
use Cwd;                      # Get the current working directory

#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $MinLength = "10";          # Min seqence length to report

#-----------------------------+
# COMMAND LINE                |
#-----------------------------+
# Get command-line arguments, or die with a usage statement
my $usage = "seq_check.pl infile \n";
my $infile = shift or die $usage;
my $infileformat = shift or die $usage;
my $CurrentDir = cwd();                    # Get the current working directory

$SeqNum = 1;

# create one SeqIO object to read in,and another to write out
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $infileformat);

 # write each entry in the input file to the new output file
while (my $inseq = $seq_in->next_seq)
{
    
    # Pad number with zeros
    $num = sprintf("%7d", $SeqNum);  
    $num =~ tr/ /0/;
    
    if ( $inseq->length < $MinLength ) {
	print "PROBLEM SEQUENCE\n";
	print "Primary ID: ".$inseq->primary_id."\n";
	print "Length: ".$inseq->length."\n";
    }
    else {
	print STDERR "Length: ".$inseq->length."\n";
    }
    
    $SeqNum++;                # Increment SeqNumber

}

exit;

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 07/11/2007
#  -Updated code format for jperl.
