#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# Name.pl                                                   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 00/00/2007                                       |
# UPDATED: 00/00/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev$                                                  |
#                                                           |
#-----------------------------------------------------------+

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
my $infile;
my $outfile;
my $inseqname;

# Booleans
my $do_gff_append = 0;
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
                    "o|outfile=s" => \$outfile,
		    "s|seqname=s" => \$inseqname,
		    # ADDITIONAL OPTIONS
		    "append"      => \$do_gff_append,
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

# SHOW HELP IF REQUIRED VARIABLES NOT PRESENT

if ( (!$infile) || (!$outfile) || (!$inseqname) ) {
    print "\a";
    print "ERROR: An input file must be specified\n" if (!$infile); 
    print "ERROR: An output file must be specified\n" if (!$outfile);
    print "ERROR: A sequence name must be specified\n" if (!$inseqname);
    print_help("full");
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
if ($do_gff_append) {
    findltr2gff ( $infile, $outfile, 1, $inseqname);
}
else {
    findltr2gff ( $infile, $outfile, 0, $inseqname);
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
	"  --seqname      # Name to use in seqname column in gff file".
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


sub findltr2gff {

    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ($findltr_in, $gff_out, $append_gff, $seqname) = @_;

    # find_ltr
    my $findltr_id;                 # Id as assigned from find_ltr.pl
    my $findltr_name;               # Full name for the find_ltr prediction
    my $ltr5_start;                 # Start of the 5' LTR
    my $ltr5_end;                   # End of the 5' LTR
    my $ltr5_len;                   # Length of the 5' LTR
    my $ltr3_start;                 # Start of the 3' LTR
    my $ltr3_end;                   # End of the 3' LTR
    my $ltr3_len;                   # Length of the 3' LTR
    my $el_len;                     # Length of the entire element
    my $mid_start;                  # Start of the LTR Mid region
    my $mid_end;                    # End of the LTR Mid region
    my $ltr_similarity;             # Percent similarity between LTRs
    my $ltr_strand;                 # Strand of the LTR

    my @in_split = ();              # Split of the infile line
    my $num_in;                     # Number of split vars in the infile

     # Initialize Counters
    my $findltr_num = 0;            # ID Number of putatitve LTR retro

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (INFILE, "<$findltr_in") ||
	die "Can not open input file:\n$findltr_in\n";

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

	my @in_split = split;
	my $num_in = @in_split;   
	
	# Load split data to vars if expected number of columns found
	if ($num_in == 10) {

	    $findltr_num++;

	    $findltr_id = $in_split[0];
	    $ltr5_start = $in_split[1];
	    $ltr5_end = $in_split[2];
	    $ltr3_start = $in_split[3];
	    $ltr3_end = $in_split[4];
	    $ltr_strand = $in_split[5];
	    $ltr5_len = $in_split[6];	    
	    $ltr3_len = $in_split[7];
	    $el_len = $in_split[8];
	    $ltr_similarity = $in_split[9];

	    $mid_start = $ltr5_end + 1;
	    $mid_end = $ltr3_start - 1;   

	    $findltr_name = $seqname."_findltr_"."".$findltr_id;

	    # 5'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"find_ltr\t".          # Source
		"exon\t".              # Features, exon for Apollo
		"$ltr5_start\t".       # Feature start
		"$ltr5_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	    # MID
	    print GFFOUT "$seqname\t". # Name of sequence
		"find_ltr\t".          # Source
		"exon\t".              # Features, exon for Apollo
		"$mid_start\t".        # Feature start
		"$mid_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)
	    
	    # 3'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"find_ltr\t".          # Source
		"exon\t".              # Features, exon for Apollo
		"$ltr3_start\t".       # Feature start
		"$ltr3_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	} # End of if num_in is 10

    } # End of while INFILE


} # End of findltr2gff

=head1 NAME

cnv_findltr2gff.pl - Convert output from find_ltr.pl to gff format.

=head1 VERSION

This documentation refers to version $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_findltr2gff.pl -i InFile.ltrpos -o OutFile.gff
    
    --infile        # Path to the ltrpos input file 
    --outfie        # Path to the gff format output file
    --seqname       # Name to use in seqname column in gff file

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

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

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
