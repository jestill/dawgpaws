#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# gff_seg.pl - Segmentation of a large gff file             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/17/2008                                       |
# UPDATED: 02/17/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a gff file that contains point or segement data    |
#  will extract segments that exceed a threshold value or   |
#  array of threshold values. Will create a gff segment     |
#  file as well as a gff parse file. The segment file       |
#  converts the vals to segments that meet the threshold    |
#  criteria while the parse file returns all points or      |
#  segments in the input file that exceed the threshold     |
#  value.                                                   |
#                                                           |
# USAGE:                                                    |
#  gff_seg.pl -i infile.gff -s outfile.gff - p parseout.gff |
#             -t [integer]
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
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
# Vars that take values from command line
my $infile;                   # Gff infile
my $thresh;                   # Threshold value for parsing, segmenting
my $outfile_seg;              # Path to the segmentation output file
my $outfile_parse;            # Path to the parsed output file
#my $outdir;                  # Base outdir for storing output when $thresh is an array

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Vars that take values from the gff file
my $seg_start;
my $seg_end;

my $in_seg = 0;               # Boolean, TRUE if in SEG

# GFF IN VALS
my $in_seq_name;  # 0 - Sequence name
my $in_source;    # 1 - Source of the feature
my $in_feat;      # 2 - feature
my $in_start;     # 3 - Start of the feature
my $in_end;       # 4 - End of the feature
my $in_score;     # 5 - Score, This is oligo count too
my $in_strand;    # 6 - Strand of the feature
my $in_frame;     # 7 - Frame of the feature
my $in_attribute; # 8 - Attribute info

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"    => \$infile,
                    "s|seg-out=s"   => \$outfile_seg,
		    "p|parse-out=s" => \$outfile_parse,
		    "t|thresh=s"    => \$thresh,
		    # ADDITIONAL OPTIONS
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
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

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if  ( (!$infile) || (!$thresh) || ( (!$outfile_seg) || ($outfile_parse)  )  ) {
    print "\a";
    print "ERROR: An infile must be specfied at the command line" 
	if (!$infile);
    print "ERROR: A threshold value must be specified at the command line"
	if (!$thresh);
}

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 

#-----------------------------+
# OPEN FILE HANDLES           |
#-----------------------------+
if ( ($outfile_seg) ) {
    open (SEGOUT, ">$outfile_seg")
	|| die "ERROR: Can not open segment out file:\n$outfile_seg\n";
}

if ( ($outfile_parse) ) {
    open (PAROUT, ">$outfile_parse")
	|| die "ERROR: Can not open parse out file:\n$outfile_parse\n";
}

open (GFFIN, "<$infile")
    || die "ERROR: Can not open gff input file:\n$infile\n";

#-----------------------------+
# PROCESS GFF FILE            |
#-----------------------------+
while (<GFFIN>) {
    
    chomp;
    my @gff_data = split(/\t/, $_);
    my $num_cols = @gff_data;
    
    # Expected number of cols is 9, otherwise
    # may want to send the output to STDERR
    # If the parse of cols is incorrect
    #print "COLS:".$num_cols."\n";
    if ($num_cols == 9) {

	# Get simple names
	$in_seq_name  = $gff_data[0];  # 0 - 
	$in_source    = $gff_data[1];  # 1 - 
	$in_feat      = $gff_data[2];  # 2 - 
	$in_start     = $gff_data[3];  # 3 - 
	$in_end       = $gff_data[4];  # 4 - 
	$in_score     = $gff_data[5];  # 5 - 
	$in_strand    = $gff_data[6];  # 6 - 
	$in_frame     = $gff_data[7];  # 7 - 
	$in_attribute = $gff_data[8];  # 8 - 

	if ( $in_score >= $thresh) {
	    
	    # The following for debug
	    #print "$in_score \:\>\: $thresh\n";

	    # PRINT PARSE DATA DIRECTLY TO OUTFILE
	    if ( ($outfile_parse) ) {
		#print STDOUT $_."\n";
		print PAROUT $_."\n";
	    }

	    #-----------------------------+
	    # SEGMENT DATA                |
	    #-----------------------------+
	    # First time we encounter a value above the threshold we are
	    # in a segment that exceeds the threshold criteria
	    if ($in_seg==0) {               # We are at the start of the segment
		$seg_start = $in_start;
		$seg_end = $in_end;
	    }
	    else {
		$seg_end = $in_end;
	    }
	    $in_seg=1;
	    
	    
	}
	else {
	    # If we just got out of a segment, print output
	    if ( $in_seg==1 ) {
		# A very simplified view of the gff file
		# print to the 
		my $gff_out_str = $in_seq_name."\t".
		    $in_source."\t".
		    $in_feat."\t".
		    $seg_start."\t".
		    $seg_end."\t".
		    $thresh."\t".
		    $in_strand."\t".
		    $in_frame."\t".
		    $in_attribute."\n";

		print SEGOUT "$gff_out_str" if ($outfile_seg);
		print STDOUT "$gff_out_str" unless ($quiet);
	    }
	    # reset in_seg to false
	    $in_seg=0;
	}

    } # End if num_cols = 9
} # End of while GFFIN



close (SEGOUT) if ($outfile_seg);
close (PAROUT) if ($outfile_parse);

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


=head1 NAME

gff_seg.pl - Segment and parse a large gff file

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS
    
  USAGE:
    gff_seg.pl -i infile.gff -s seg_out.gff -p par_out.gff -t integer

    -i,--infile         # Path to the input file
    -s,--seg-out        # Path to the segmented output file
    -p,--parse-out      # Path to the parsed output file
    -t,--thresh         # Threshold value

=head1 DESCRIPTION

Given a gff file that contains point or segement data
will extract segments that exceed a threshold value or
array of threshold values. Will create a gff segment
file as well as a gff parse file. The segment file
converts the vals to segments that meet the threshold
criteria while the parse file returns all points or
segments in the input file that exceed the threshold
value.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -s,--seg-out

Gff file of the segmented data

=item -p,--parse-out

Gff file of the parsed data. This will include all rows in the original
file that exceed the threshold value.

=item -t,--thresh

Threshold value. A single integer that represents the threshold value.
The output will return all features that are greater then or equal
to the threshold value.

=back

=head1 Additional Options

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

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 02/17/2008

UPDATED: 02/17/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/17/2008
# - Starting program. Treating thresh as a single value
#   later updates will need to accept an array here.
# - ASSUMES CONTIGUOUS FEATURES
# - Overlap in segments okay
# - Two output files will be created:
#     - parse_out - All features exceeding the threshold value
#     - seg_out - All segments exceeding the threshold value
