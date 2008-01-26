#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq2dir.pl - Split sequence file to dir of fasta files    | 
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/22/2004                                       |
# UPDATED: 01/26/2008                                       |
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

    my $seq_id = $inseq->primary_id();
    
    # Pad number with zeros so that the total length of the 
    # string is 7 characaters (ie. 0000012)
    $num = sprintf("%7d", $SeqNum);
    $num=~ tr/ /0/;
    my $out_file_path = $OutputDir.$seq_id.".fasta";

    my $seq_out = Bio::SeqIO->new('-file' => ">$out_file_path",
				  '-format' => fasta);
    print STDERR "Processing: $seq_id\n";
    #print $outfile.$num." \n";
    
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
# 01/26/2008
#  - Changed default outfile name to the header name for
#    the sequence file from the fasta record
