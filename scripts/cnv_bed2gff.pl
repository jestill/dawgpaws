#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_bed2gff.pl - Convert BED file to GFF format.          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/01/2011                                       |
# UPDATED: 08/01/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert BED file to GFF format.                          |
#                                                           |
# USAGE:                                                    |
#  cnv_bed2gff.pl -i infile.bed -o outfile.gff              |
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
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

my $thresh_val = 0;               # Threshold value
my $seqname;                      # Sequence name to search for
my $do_append=0;                  # Append result to existing GFF file
my $prog = "TopHat";              # The name of the program
my $param;                        # The parameter set, db name used
my $delim = ":";                  # Delimiting character to use in GFF outfile

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "delim=s"     => \$delim,
		    "program=s"   => \$prog,
		    "p|param=s"   => \$param,
		    "n|s|seqname|name=s"  => \$seqname,
		    "gff-ver=s"   => \$gff_ver,
		    "t|thresh=i"  => \$thresh_val,
		    "append"      => \$do_append,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);


#-----------------------------+
# STANDARDIZE GFF VERSION     |
#-----------------------------+
unless ($gff_ver =~ "GFF3" || 
	$gff_ver =~ "GFF2") {
    # Attempt to standardize GFF format names
    if ($gff_ver =~ "3") {
	$gff_ver = "GFF3";
    }
    elsif ($gff_ver =~ "2") {
	$gff_ver = "GFF2";
    }
    else {
	print "\a";
	die "The gff-version \'$gff_ver\' is not recognized\n".
	    "The options GFF2 or GFF3 are supported\n";
    }
}

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
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
    print "\ncnv_bed2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
bed2gff ($infile,
	 $prog,
	 $param,
	 $delim,
	 $outfile,
	 "GFF2",
	 $seqname,
	 $thresh_val,
	 $do_append);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub bed2gff {

    my ($bed_in_file,
	$prog,
	$param,
	$delim,
	$gff_out,
	$gff_format,
	$seqname,
	$thresh,
	$do_append)  =  @_;
    
    
    # Append parameter name to program name
    # ie TOPHAT:Flower or TOPHAT:Stem
    if ($param) {
	$prog = $prog.$delim.$param;
    }

    #-----------------------------+
    # OPEN BED FILE               |
    #-----------------------------+
    # Default is to expect intput from STDIN if infile path not given
    if ($bed_in_file) {
	open ( BED_IN, "<".$bed_in_file ) ||
	    die "Could not open the BED infile:\n $bed_in_file\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open ( BED_IN, "<&STDIN" ) ||
	    die "Could not open STDIN for input\n";
    }
    
    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gff_out) {
	if ($do_append) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
	else {
	    open (GFFOUT,">$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	    if ($gff_ver =~ "GFF3") {
		print GFFOUT "##gff-version 3\n";
	    }
	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
	if ($gff_ver =~ "GFF3") {
	    print GFFOUT "##gff-version 3\n";
	}
    }
    
    
    while (<BED_IN>) {
	chomp;

	# GET THE PARTS OF THE BED INPUT FILE
	my @bed_parts = split;
	my $bed_chrom = $bed_parts[0];
	# Bed starts chromosome at 0 so we need to add one
	my $bed_start = $bed_parts[1] + 1;
	my $bed_end = $bed_parts[2];
	my $bed_name = $bed_parts[3];
	my $bed_score = $bed_parts[4];
	my $bed_strand = $bed_parts[5];
	my $bed_thick_start = $bed_parts[6];
	my $bed_thick_end = $bed_parts[7];
	my $bed_reserved = $bed_parts[8];
	my $bed_block_count = $bed_parts[9];
	my $bed_block_sizes = $bed_parts[10];
	my $bed_block_starts = $bed_parts[11];
	my $bed_exp_count = $bed_parts[12];
	# The following is a comma separated list
	my $bed_exp_ids = $bed_parts[13];
	# The following is a comma separated list
	my $bed_exp_scores = $bed_parts[14];


	#-----------------------------+
	# SKIP VALUES BELOW THRESHOLD |
	#-----------------------------+
	if ( $thresh ) {
	    if ($bed_score < $thresh) {
		next;
	    }
	}

	#-----------------------------+
	# SKIP VALUES NOT IN THE      |
	# CONTIG WE ARE LOOKING FOR   |
 	#-----------------------------+
	if ($seqname) {
	    unless ($seqname =~ $bed_chrom) {
		next;
	    }

	}


	# Start to print GFF output
	my $gff_seqid = $bed_chrom;

	#-----------------------------+
	# SET THE GFF ATTRIBUTE       |
	#-----------------------------+
	# Need to be able to set different attirbute for GFF2 or GFF3 format
	my $gff_attribute;
	if ($gff_ver =~ "GFF2") {
	    $gff_attribute = $bed_name;
	}
	elsif ($gff_ver =~ "GFF3") {
	    $gff_attribute = "ID=".$bed_name;
	}
	else {
	    print STDERR "WARNING: GFF VERSION ".$gff_ver.
		" NOT RECOGNIZED!!\n";
	}

	my $gff_str = $gff_seqid."\t".   # Contig name
	    $prog."\t".                  # Source (Program Name)
	    "match"."\t".                # Feature type ie match
	    $bed_start."\t".             # Feature start
	    $bed_end."\t".               # Feature end
	    $bed_score."\t".             # Score
	    $bed_strand."\t".            # Strand
	    ".\t".                       # Phase
	    $gff_attribute;


	# Print GFF STRING if meets criteria
	print GFFOUT $gff_str."\n";
	    
    }

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

1;
__END__

=head1 NAME

cnv_bed2gff.pl - Convert BED file to GFF format. 

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_bed2gff.pl -i infile.bed -o outfile.gff

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

Convert BED format files to GFF format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 OPTIONS

=over 2

=item --program

The name of the program to use for the third column of the GFF
output file. By default this is set to TopHat.

=item --param

The name of the paramter set to tag the third column of the gff outpt 
file with. Be default this is blank. This parameter name will be
added to the program name in the third column of the GFF output file. 

=item --delim

The delimiting character to use to separate the program from the 
paramter name. For Apollo it is convienient to use the : character but 
this causes errors in gbrowse. This allows differenct characters such as'_' 
to be used for GFF output.

=item -n,--name

The name of the sequence contig to return results for. For example, in large
BED files including results from multiple contigs we may just want to
return values for a single contig.

=item --gff-ver

The GFF version to produce in output. Valid values are GFF2 or GFF3.

=item -t,--thresh

The threshold value from the BED file to return values for. All values equal
to or greater than this threshold value will be returned in the GFF output
file. This for example allows for a filter to to ignore TopHat junctions
that are only returned by a single match.

=item --append

If an file exists at the path name by --outfile, this will allow for 
any results to be appended to the existing outfile.

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

An example use is given the junctions.bed output from the TopHat
program to extract values from a single scaffold that occur more
than one time. The GFF output file with use the _ delimiting
character to facilitate dispaly in the Gbrowse program.

  cnv_bed2gff.pl -i junctions.bed -o scaffold00001_junctions.gff
                 --thresh 1 -n scaffold00001 --delim _
                 --gff-ver GFF2

Output from this file would be similar to:

  scaffold00001	TopHat_min3	match	29346	43483	3	+	.	JUNC00084836
  scaffold00001	TopHat_min3	match	85197	85461	125	-	.	JUNC00084837
  scaffold00001	TopHat_min3	match	85403	85655	171	-	.	JUNC00084839
  scaffold00001	TopHat_min3	match	85677	105989	207	-	.	JUNC00084840
  scaffold00001	TopHat_min3	match	106075	107480	231	-	.	JUNC00084843
  scaffold00001	TopHat_min3	match	107488	107799	281	-	.	JUNC00084844
  scaffold00001	TopHat_min3	match	107493	107775	6	-	.	JUNC00084845
  scaffold00001	TopHat_min3	match	107728	113574	168	-	.	JUNC00084847
  scaffold00001	TopHat_min3	match	107764	113581	3	-	.	JUNC00084848
  scaffold00001	TopHat_min3	match	107795	113558	5	-	.	JUNC00084851

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

The following variables may be set in the user environment.

=over 2

=item DP_GFF

The deafult value to use for the GFF output format. Valid values 
are GFF2 or GFF3.

=back

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=head1 REFERENCE

Your use of the DAWGPAWS programs should reference the following manuscript:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

Additional information on the BED format are available at:
http://genome.ucsc.edu/FAQ/FAQformat.html#format1

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/01/2011

UPDATED: 08/01/2011

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 08/01/2011
# - Started the program for basic parsing of BED format
#   files to the GFF format for the purpose of parsing
#   TopHat output for use with Gbrose or for use as 
#   input into the Augustus genome annotation program.
