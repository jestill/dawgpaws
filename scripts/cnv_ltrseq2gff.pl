#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# ltr_seq2gff.pl - Convert LTR_seq output to gff format.    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/03/2007                                       |
# UPDATED: 09/03/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts LTR_seq output to gff file format.              |
#                                                           |
#-----------------------------------------------------------+
# TODO: Extract only the unique LTR Predictions 

=head1 NAME

cnv_ltrseq2gff.pl - Convert LTR_seq output to gff format.

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

 Usage:
  Name.pl -i InFile -o OutFile

=head1 DESCRIPTION

Converts program output from LTR_seq to 

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 OPTIONS

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

B<LTR_seq Program:>

This program requires output from LTR_seq. LTR_seq is available upon email
request from the author, Ananth Kalyanaraman: 
http://www.eecs.wsu.edu/~ananth/contact.htm.

Also see the original publication for the LTR_par program :

I<A. Kalyanaraman, S. Aluru. 2006. Efficient algorithms and software for 
detection of full-length LTR retrotransposons. Journal of 
Bioinformatics and Computational Biology (JBCB), 4(2):197-216.>

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Path of the LTR_seq file to convert
my $outfile;                   # Path to the gff format file created
my $inseqname;                 # Name of the sequence analyzed

# Optional
my $infasta;                   # Path to the FASTA file analyzed

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_append = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required options
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    "s|seqname=s"   => \$inseqname,
		    # Additional Optons
		    "f|fastafile=s" => \$infasta,
		    "append"        => \$do_gff_append,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # Additional information
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}

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

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

if ($do_gff_append) {
    ltrseq2gff ( $infile, $outfile, 1, $inseqname);
} 
else {
    ltrseq2gff ( $infile, $outfile, 0, $inseqname);
}

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

sub ltrseq2gff {

    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ($ltrseq_in, $gff_out, $append_gff, $seqname) = @_;

    #-----------------------------+
    # LTRSEQ VARS                 |
    #-----------------------------+
    my $ltrseq_infile;              # Input file that was analyzed by LTR_seq
    my $ltrseq_numseqs;             # Number of seqs analyzed by LTR_seq
    my $ltrseq_msrt;                # Max Score Ration Threshold
    my $ltrseq_ltrmin;              # LTRmin
    my $ltrseq_ltr_minexact;        # LTR Min Exact Match
    my $ltrseq_dmin;                # Dmin
    my $ltrseq_dmax;                # Dmax

    # Postion of the LTR Retrotransposon features
    my $ltrseq_name;                # Unique name assigned to the prediction
    my $ltrspan_start;              # Start of the entire LTR Retrotransposon
    my $ltrspan_end;                # End of the entire LTR Retrotransposon
    my $ltrspan_len;                # Length of the LTR span
    my $ltr5_start;                 # Start of the 5' LTR
    my $ltr5_end;                   # End of the 5' LTR
    my $ltr3_start;                 # Start of the 3' LTR
    my $ltr3_end;                   # End of the 3' LTR
    my $ltr_len;                    # Length of the LTRs
    my $mid_start;                  # Start of the LTR Mid region
    my $mid_end;                    # End of the LTR Mid region
    my $ltr_diff;                   # Percent Difference in the LTRs
    my $ltr_conf;                   # Confidence score for the prediction
    my $ltr_tsr;                    # Target site rep

    my @in_split = ();              # Split of the infile line
    my $num_in;                     # Number of split vars in the infile

    # Initialize Counters
    my $ltrseq_num = 0;             # ID Number of putatitve LTR retro (Incremented Number)

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (INFILE, "<$ltrseq_in") ||
	die "Can not open input file:\n$ltrseq_in\n";

    if ($append_gff) {
	open (GFFOUT, ">>$gff_out") ||
	    die "Could not open output file for appending\n$gff_out\n";
    }
    else {
	open (GFFOUT, ">$gff_out") ||
	    die "Could not open output file for output\n$gff_out\n";
    } # End of if append_gff

    #-----------------------------+
    # PROCESS INFILE              |
    #-----------------------------+
    while (<INFILE>) {
	chomp;

	# Print the INFILE
	#print "$_\n";

	if (m/^Report:/) {

	    # If this is a Report line, accepted or rejected can be
	    # determined by counting the line number
	    # 15 is a rejected line
	    # 19 is an accepted line 
	    my @in_split = split;
	    my $num_in = @in_split;   

	    #print "Report Line: $num_in parts\n";
	    #print "\t$_\n";


	    if ($num_in == 15) {
		#print "\tREJECTED\n";
	    }
	    #-----------------------------+
	    # ACCEPTED LTR RETRO Parts    | 
	    #-----------------------------+
	    elsif ($num_in == 19) {
		#print "\tACCEPTED\n";
		print "\n".$_."\n";
		$ltrseq_num++;              # Increment the count

		$ltr5_start = $in_split[4] || "NULL";
		$ltr5_end = $in_split[5] || "NULL";
		$ltr3_start = $in_split[7] || "NULL";
		$ltr3_end = $in_split[8] || "NULL";
		$mid_start = $ltr5_end + 1;
		$mid_end = $ltr3_start - 1;
		$ltrspan_start = $ltr5_start;
		$ltrspan_end = $ltr3_end;
		$ltrspan_len = $in_split[12];
		$ltr_len = $in_split[10];
		$ltr_diff = $in_split[14];
		$ltr_conf = $in_split[16];
		$ltr_tsr = $in_split[18];

		$ltrseq_name = $seqname."_ltrseq_".$ltrseq_num;

		if ($ltr_tsr =~ m/\#(.*)\#/) {
		    $ltr_tsr = $1;
		}

		# SHOW INFO IF VERBOSE
		
		if ($verbose) {
		    print STDERR "$ltrseq_name\n";
		    print STDERR "\tSTART:    $ltrspan_start\n";
		    print STDERR "\tEND:      $ltrspan_end\n";
		    print STDERR "\tLTRLEN:   $ltr_len\n";
		    print STDERR "\tSPAN_LEN: $ltrspan_len\n";
		    print STDERR "\tLTSR:     $ltr_tsr\n";
		    print STDERR "\tLTRCON:   $ltr_conf\n";
		    print STDERR "\tLTRDIF:   $ltr_diff\n";
		} # End of if verbose

		#-----------------------------+
		# PRINT TO GFF OUTPUT FILE    |
		#-----------------------------+
		# May want to allow choice for showing
		# Details or not

#		# LTR RETRO PREDICTION SPAN
#		print GFFOUT "$seqname\t".     # Name of sequence
#		    "LTR_seq:span\t".          # Source name
#		    "exon\t".                  # Feature, exon for Apollo
#		    "$ltrspan_start\t".        # Start of the ltr span
#		    "$ltrspan_end\t".          # End of the ltr span
#		    "$ltr_conf\t".             # Score, LTR Confidence Score
#		    ".\t".                     # Strand
#		    ".\t".                     # Frame
#		    "$ltrseq_name"."_span\n";  # Features (Name)
		
		
		# 5' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "LTR_seq\t".               # Source
		    "exon\t".                  # Feature, exon for Apollo
		    "$ltr5_start\t".           # Start of the 5'ltr
		    "$ltr5_end\t".             # End of the 5' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

		# MID
		print GFFOUT "$seqname\t".     # Name of sequence
		    "LTR_seq\t".               # Source
		    "exon\t".                  # Feature, exon for Apollo
		    "$mid_start\t".            # Start of the 3'ltr
		    "$mid_end\t".              # End of the 3' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

		# 3' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "LTR_seq\t".               # Source
		    "exon\t".                  # Feature, exon for Apollo
		    "$ltr3_start\t".           # Start of the 3'ltr
		    "$ltr3_end\t".             # End of the 3' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

	    }


	    # 15 Parts is Rejected

	}
	# HEADER INFORMATION
	elsif (m/^Info:/) {
	    
	}
	# The input string that was passed
	elsif (m/^Input:/) {

	}

    } # End of while INFILE

}

=head1 HISTORY

STARTED: 09/03/2007

UPDATED: 09/03/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/03/2007
