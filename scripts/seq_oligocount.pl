#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq_oligocount.pl - Count oligos of size k
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/11/2007
# UPDATED: 10/12/2007
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;                # Seq IO used to to work with the input file
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Path for the input query sequence
my $outdir;                    # Path for the output
my $index;                     # Path to the sequence index file
my $kmer_len = 20;             # Oligomer length, default is 20
my $seq_name;                  # Sequence name used in the output

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
		    "d|db=s"      => \$index,
                    "o|outdir=s"  => \$outdir,
		    "n|name=s"    => \$seq_name,
		    # ADDITIONAL OPTIONS
		    "k|kmer=s"    => \$kmer_len,   
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
#if ($show_usage) {
#    print_help("");
#}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\n$0:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

# Throw error if required options not present
if ( (!$infile) || (!$outdir) || (!$index) ) {

    print "\a";
    print "ERROR: Input file path required" if !$infile;
    print "ERROR: Output directory required" if !$outdir;
    print "ERROR: Index file path required" if !$index;
    print_help("");

}

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 

seq_kmer_count ($infile,$outdir,$index,$kmer_len, $seq_name);

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}


sub seq_kmer_count {

    my ($fasta_in, $outdir, $vmatch_index, $k, $seq_name) = @_;
    
    # Array of threshold values
    #my @thresh = ("200");
    my $thresh = 50;

    my $in_seq_num = 0;
    my $inseq = Bio::SeqIO->new( -file => "<$fasta_in",
				 -format => 'fasta');

    # Counts of the number of occurences of each oligo
    my @counts = ();    
    my @vdat;   # temp array to hold the vmatch data split by tab
    my $pos = 0;

    my $start_pos;
    my $end_pos;
    my $i;         # Array index value

    # FILE PATHS
    # These will need to be modified later to more useful names
    my $temp_fasta = $outdir."split_seq_".$k."_mer.fasta";
    my $vmatch_out = $outdir."vmatch_out.txt";
    my $gff_count_out = $outdir."vmatch_out.gff";

    while (my $seq = $inseq->next_seq) {

	$in_seq_num++;
	if ($in_seq_num == 2) {
	    print "\a";
	    die "Input file should be a single sequence record\n";
	}

	# Calculate base cooridate data
	my $seq_len = $seq->length();
	my $max_start = $seq->length() - $k;
	
	# Print some summary data
	print STDERR "\n==============================\n" if $verbose;
	print STDERR "SEQ LEN: $seq_len\n" if $verbose;
	print STDERR "MAX START: $max_start\n" if $verbose;
	print STDERR "==============================\n" if $verbose;
	
	# CREATE FASTA FILE OF ALL K LENGTH OLIGOS
	# IN THE INPUT SEQUENCE
	print STDERR "Creating oligo fasta file\n" if $verbose;
	open (FASTAOUT, ">$temp_fasta") ||
	    die "Can not open temp fasta file:\n $temp_fasta\n";

	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;

	    my $oligo = $seq->subseq($start_pos, $end_pos);

	    # Set counts array to zero
	    $counts[$i] = 0;
	    
	    print FASTAOUT ">$start_pos\n";
	    print FASTAOUT "$oligo\n";
	    
	}

	close (FASTAOUT);

	# QUERY OLIGO FASTA FILE AGAINST VMATCH INDEX
	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
	    " $vmatch_index > $vmatch_out";
#	# OUTPUT TO STDOUT FOR DEBUG
#	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
#	    " $vmatch_index";
	print STDERR "\nVmatch cmd:\n$vmatch_cmd\n" if $verbose;
	system ($vmatch_cmd);

	# PARSE VMATCH OUTPUT FILE
	# increment
	unless (-e $vmatch_out) {
	    print STDERR "Can not file the expected vmatch output file\n";
	    die;
	}

	# PARSE THE VMATCH OUTPUT FILE AND INCREMENT
	# THE COUNTS ARRAY
	open (VMATCH, "<$vmatch_out") ||
	    die "Can not open vmatch output file:\n$vmatch_out\n";

	# Count oligos hits in index file
	print STDERR "\nCounting oligos ...\n" if $verbose;
	while (<VMATCH>) {
	    # ignore comment lines
	    unless (m/^\#/) {
		chomp;
		@vdat = split;
		

		my $len_vdat = @vdat;
		# print the following for debu
		#print "VDAT LEN: $len_vdat\n";
		#print $vdat[5]."\n";
		
		# Get the seq file in the new sub se
		# Counts index starts at zero
		# It may be possible to increment without needing
		# to add one
		#$counts[$vdat[5]]++;
		$counts[$vdat[5]] = $counts[$vdat[5]] + 1;

	    }
	}
	close (VMATCH);


	#///////////////////////////////////
	# NEED SEGMENTATION STEP HERE
	# THIS WILL JOIN INDIVIDUAL PIPS
	# FOR A GIVEN THRESHOLD COVERAGE
	# This could also be done with a 
	# separate script operating on the
	# gff output from this program.
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	# Output oligo counts to gff file
	
	# This is the pip file, a single pip with depth coverage data
	# for eack oligo
	open (GFFCOUNT, ">$gff_count_out") ||
	    die "Can not open gff out file:\n$gff_count_out\n";
	

	print "\nCreating gff output files ...\n";
	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;
	    
	    ## print output to stdout
	    #print GFFCOUNT "$start_pos\t".$counts[$i]."\n";
	    
	    # The count will be placed in the score position
	    my $seq_str = $seq->subseq($start_pos, $end_pos);
	    print GFFCOUNT "$seq_name\t".  # Ref sequence
		"vmatch\t".                # Source
 		"wheat_count\t".           # Type
		"$start_pos\t".            # Start
		"$end_pos\t".              # End
		$counts[$i]."\t".          # Score
		".\t".                     # Strand
		".\t".                     # Phase
		"Vmatch ".$k."mer\n";    # Group


	    # PRINT OUT SEQS EXCEEDING THRESHOLD VALUES
	    if ($counts[$i] > $thresh) {
		my $thresh_seq = $seq->subseq($start_pos, $end_pos);
		print "$start_pos\t".$counts[$i]."\t".
		    "$thresh_seq\n";
	    }

	}
	
	close (GFFCOUNT);

    } # End of while seq object

    # May want to make a single fasta file of all oligos for the
    # qry fasta_sequence and parse the results from vmatch from that
    

}

=head1 NAME

seq_oligocount.pl - Count oligos from an input sequence

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head 2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

=over

=over 2

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

=head2 Perl Modules

This program requires the following Perl modules

B<Bio:SeqIO>

=head2 Software

B<vmatch>

This program requires the vmatch series of programs.
H<http://www.vmatch.de>

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
