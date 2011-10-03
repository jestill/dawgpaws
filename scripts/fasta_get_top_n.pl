#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_get_top_n.pl - Get the first n seqs of a fasta file |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/17/2011                                       |
# UPDATED: 02/17/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Goal is to fetch the top n seqences of a fasta file.     | 
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
use Getopt::Long;             # Get options from the command line

#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $min_length = int("10");    # Min seqence length to report
my $num_seqs = int('5');       # The path to the intput file
my $infile;                    # The path to the intput file
my $outfile;                   # The path to the output file
my $infile_format = 'fasta';   # The format of the input sequence file
my $outfile_format = $infile_format ;

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $is_gap = 0;                # Current character matche the gap char
my $prev_gap = 0;              # Previous character matches the gap char

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
                    "n|num-seqs=i"  => \$num_seqs, 
                    # Additional options
                    "f|in-format=s" => \$infile_format,
                    "out-format-s"  => \$outfile_format,
		    "l|length=i"    => \$min_len,
		    # Booleans
		    "verbose"       => \$verbose,
		    "test"          => \$test,
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
                    "h|help"        => \$show_help, );

# Get command-line arguments, or die with a usage statement
my $usage = "fasta_get_top_n.pl -i infile -o outfile -l 100000 -f fasta";


# create one SeqIO object to read in,and another to write out
my $seq_in;

# If no input file path is given expects input from STDIN
if ($infile) {
    $seq_in = Bio::SeqIO->new(-file   => "<$infile",
			      -format => $infile_format);
}
else {
    $seq_in = Bio::SeqIO->new(-fh     => \*STDIN,
			      -format => $infile_format);
}

# If not outfile path, write output to STDOUT
my $seq_out;
if ($outfile) {
    $seq_out = Bio::SeqIO->new(-file   => ">$outfile",
			       -format => $outfile_format);
}
else {
    $seq_out = Bio::SeqIO->new(-fh     => \*STDOUT,
			      -format  => $outfile_format);
}


 # write each entry in the input file to the new output file
my $seq_count = 0;
while (my $inseq = $seq_in->next_seq) {
 
#    if ($min_len) {
#	if ( $inseq->length < $min_length ) {
#	    next;
#	}
#    }
    
    $seq_out->write_seq($inseq);
    $seq_count++;

    last if $seq_count == $num_seqs;
    
}

print STDERR "$seq_count sequences were processed\n";

exit;

__END__

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 02/17/2011
# Program started

