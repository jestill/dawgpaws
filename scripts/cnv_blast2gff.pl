#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2gff.pl - Convert BLAST output to gff            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 01/31/2009                                       |
#                                                           |  
# DESCRIPTION:                                              | 
#  Convert blast output to a Apollo compatible gff file.    |
#                                                           |
# USAGE:                                                    |
#  cnv_blast2gff.pl -i infile.bln -o blast_out.gff          |
#                                                           |
# VERSION:                                                  |
#  $Rev$                                              |
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
use Bio::SearchIO;             # Parse BLAST output
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $max_e;
my $min_len;

# Set variable scope
#my $seqname = "seq";
my $qry_name;
my $infile;                     # Default to STDIN if not given
my $outfile;                    # Default to STDOUT if not given
my $blast_program = "blast";
my $blast_alignment = 0;        # Alternatives include 8 and 9
my $param;

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    "n|name=s"     => \$qry_name,
		    # OPTIONS
		    "p|program=s"  => \$blast_program,
		    "m|align=s"    => \$blast_alignment,
		    "d|database=s" => \$param,
		    "e|maxe=s"     => \$max_e,
		    "l|minlen=s"   => \$min_len,
		    # Booleans
		    "verbose"      => \$verbose,
		    "append"       => \$do_append,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);
 
#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( $show_usage ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\ncnv_blast2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
blast2gff ($infile, $outfile, $do_append, $qry_name, $blast_alignment,
	   $param, $blast_program);

exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    # align   - the alignment type of the blast output
    # suffix  - the suffix to add to type of blast
    # prog    - blast program (blastn, blastx, wublast etc)
    # gff second column is built as progparam
    my ($blastin, $gffout, $append, $seqname, $align, $suffix, $prog) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    # Recode from number for tab delimited blast
    if ($align == 8 || $align ==9) {
	$align = "tab";
    }
    
    #-----------------------------+
    # BLAST INTPUT OBJECT         |
    #-----------------------------+
    # This allows for input from NCBI-BLAST, and WUBLAST in either 
    # default or tab delimited output as well as 
    # INPUT FROM FILE PATH
    my $blast_report;
    if ($blastin) {
	
	# TAB ALIGNED BLAST OUTPUT
	if ($align =~ "tab" ) {
	    $blast_report = new Bio::SearchIO::blasttable ( 
		'-format' => 'blasttable',
		'-file'   => $blastin )
		|| die "Could not open BLAST input file:\n$blastin.\n";
	}
	else {
	    $blast_report = new Bio::SearchIO ( '-format' => 'blast',
						'-file'   => $blastin,
						'-signif' => $max_e,
						'-min_query_len' => $min_len) 
		|| die "Could not open BLAST input file:\n$blastin.\n";
	} # Default to default blast alignment
    }

    # INPUT FROM STDIN
    else {
	if ($align =~ "tab") {
	    $blast_report = new Bio::SearchIO::blasttable ( 
		'-format' => 'blasttable',
		'-fh'     => \*STDIN )
		|| die "Could not open STDIN for blast table input\n";
	}
	else {
	    print STDERR "Expecting input from STDIN\n";
	    $blast_report = new Bio::SearchIO ( '-format' => 'blast',
						'-fh'     => \*STDIN,
						'-signif' => $max_e,
						'-min_query_len' => $min_len) 
		|| die "Could not open STDIN for blast input\n";
	}
    }
    
    #-----------------------------+
    # FILEHANDLE FOR OUTPUT       |
    #-----------------------------+
    if ($gffout) {
	if ($append) {
	    open (GFFOUT, ">>$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
	else {
	    open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }
    
    while (my $blast_result = $blast_report->next_result())
    {

	# Not using the following for the moment
	# 01/30/2009
	$blastprog = $blast_result->algorithm;
	#$dbname = $blast_result->database_name();

    	while (my $blast_hit = $blast_result->next_hit())
	{

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		my $hitname = $blast_hit->name();

		$strand = $blast_hsp->strand('query');
		
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    $strand = "+";
		}
		
		#-----------------------------+
		# SET START AND END           |
		#-----------------------------+
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();
		    $end = $blast_hsp->end();
		}
		else {
		    #
		    $start = $blast_hsp->end();
		    $end = $blast_hsp->start();
		}

		#-----------------------------+
		# GET ALGORITHM IF UNKNOWN    |
		#-----------------------------+
		# This attempts to get the algorithm from the
		# blast report, otherwise the generic 'blast' is used
		unless ($prog) {
		    if ($blast_hit->algorithm) {
			$prog = $blast_hit->algorithm;
		    }
		    else {
			$prog = "blast";
		    }
		}

		#-----------------------------+
		# APPEND DB IDENTIFIER        |
		#-----------------------------+
		# Providing a suffix at the command line will override
		# the name provided in the blast report
		my $source = $prog;
		if ($suffix) {
		    $source = $source.":".$suffix;
		}
		elsif ($blast_result->database_name()) {
		    $source = $source.":".$blast_result->database_name();
		}


		#-----------------------------+
		# GET SEQNAME IF NOT PROVIDED |
		#-----------------------------+
		# First attempt to extract from blast report
		# otherwise use 'seq' as the indentifier.
		# For multiple query sequence, this will use
		# the same identifier for all
		unless ($seqname) {
		    if ($blast_result->query_name) {
			$seqname = $blast_result->query_name();
		    }
		    else {
			$seqname = "seq";
		    }
		}
		
		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		print GFFOUT "$seqname\t".                   # Seqname
		    "$source\t".                             # Source
		    "exon\t".                                # Feature type name
		    "$start\t".                              # Start
		    "$end\t".                                # End
		    $blast_hsp->score()."\t".                # Score
		    "$strand\t".                             # Strand
		    ".\t".                                   # Frame
		    "$hitname\n";                            # Feature name

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
    close GFFOUT;
    
}

sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

=head1 NAME

cnv_blast2gff.pl - Convert blast output to GFF

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_blast2gff.pl -i blast_report.bln -o blast_out.gff -n seq_id

=head2 Required Arguments

    -i      # Path to the input file
            # If not specified the program will expect input from STDIN
    -o      # Path to the output file
            # If not specified the program will write to STDOUT
    -n      # Name of the query sequence
            # If a name is not specified, the name will be
            # extracted from the blast report or use "seq"

=head1 DESCRIPTION

This program will translate a blast report for a single query sequence into
the GFF format. Since this uses the general bioperl BLAST parser code,
this should also be able to parse output from BLAT or wublast. This code
works best for converting a blast report of a single query sequence against
a single database.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not specified, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file. If an output file is not specified, the program
will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item -n,--name

Identifier of the sequence that has been used as the query sequence
in the blast report. This will be used as the first column in the 
gff output file.

=item -p,--program

The blast program used. This will be used to identify the source program in the
second column of the GFF output file. Example of valid values include
blastn, blastx, or wublast.

=item -d,--database

The name for the database that was blasted against. If provided, this
will be appended to the program variable in the second colum of the 
GFF output file.

=item -m,--align

The alignment format use in the BLAST report to be parsed.
The program will assume that you are using the default alginment format for
blast. Otherwise, you can specify 'tab' or '8' or '9' for tab delimited blast.

=item -e,--maxe

The maximum e value threshold to accept.

=item -l,--min-len

The minimum length to accept.

=item --verbose

Run the program with maximum reporting of error and status messages.

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

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or any
variables set in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

This program requires output from a BLAST program. Since this program makes
use of the BioPerl blast parser, it should be possible to convert
local alignment results from any of the following programs:

=over 2

=item * BLAST

=item * PSIBLAST

=item * PSITBLASTN

=item * RPSBLAST

=item * WUBLAST

=item * bl2seq

=item * WU-BLAST

=item * BLASTZ

=item * BLAT

=item * Paracel 

=item * BTK

=back

=head2 Required Perl Modules

The following perl modules are required for this program:

=over 2

=item * Bio::SearchIO

The SearchIO module is part of the bioperl module.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited BLAST testing

This program has only been tested with BLAST output from the NCBI-BLAST package
using standard and tab delimited output. If you find that there are limiations
to this program that limit your use, please email the author or file a 
bug report on the Sourceforge website:
http://sourceforge.net/tracker/?group_id=204962

=back

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 01/31/2009

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started
# -Main body of program written
# 01/30/2009
# -Updated POD documentation
# -Moved POD documentation to end of the program
# -Updated to new print_help subfunction that extracts
#  help from POD
# -Added default to print to STDOUT if not outfile
#  path given at command line
# -Added default input from STDIN in no infile path
#  given at command line
# -Added package name DAWGPAWS
# -Added options to pass database name, query sequence id,
#  and program used. Otherwise these are attempted to
#  be extracted from the command line
# 01/31/2009
# -Fixed print_help to extract help from POD
