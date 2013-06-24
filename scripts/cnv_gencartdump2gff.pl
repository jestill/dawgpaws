#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gencartdump2gff.pl - Gencarout feature dump to gff3   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/08/2013                                       |
# UPDATED: 05/08/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts a feature dump of non-overlapping TEs from the  |
#  GenCART schema to a GFF3 format file.                    |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
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
my ($VERSION) = q$Rev$ =~ /(\d+)/ || "pre-release";
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

my $gff_source;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "source=s"    => \$gff_source,
		    "gff-ver=s"   => \$gff_ver,
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
# PRINT REQUESTED HELP        |
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

if ($infile) {
    open (INFILE, "<$infile") ||
	die "Can not open input file $infile";
}
else {
    open (INFILE, "<&STDIN") ||
	die "Can not open STDIN for input"
}

if ($outfile) {
    open (GFFOUT, ">$outfile") ||
	die "Can not open outfile $outfile";
}
else {
    open (GFFOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output";
}

print GFFOUT "##gff-version 3\n";

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
my $line_num;
while ( <INFILE> ) {
    chomp;
    $line_num++;
    print STDERR $_."\n"
	if $verbose;

    my @inparts = split ( /\t/, $_);
    my $scaffold_num = sprintf( "%05d", $inparts[1] );
    my $scaffold = "AmTr_v1.0_scaffold".$scaffold_num;
    my $exemplar = $inparts[11];
    my @exemplar_parts = split (/\_/, $exemplar);
    my $type = $exemplar_parts[0];
    my $start = $inparts[4];
    my $end = $inparts[6];
    my $score = $inparts[7];
    my $strand = $inparts[9];
    my $phase = ".";
    my $attribute = "ID=".$gff_source."_".$exemplar."_".
	$scaffold."_".$start."_".$end.
	";Name=".$exemplar;
    
    print GFFOUT $scaffold."\t".
	$gff_source."\t".
	$type."\t".
	$start."\t".
	$end."\t".
	$score."\t".
	$strand."\t".
	$phase."\t".
	$attribute.
	"\n";
	
}

close (INFILE);
close (GFFOUT);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    cnv_gencartdump2gff.pl -i in_sql_dump.txt -o outfile.gff3

=head2 Required Arguments

    --infile        # Path to the input file from SQL dump
    --outfie        # Path to the output file

=head1 DESCRIPTION

Converts a feature dump of non-overlapping TEs from the
GenCart schema to a GFF3 format file suitable for upload
to another genome feature database.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file which is a text file resulting from a mysqldump
from the Gencart database.

=item -o,--outfile

Path of the output file, which will be a GFF3 format file.

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

A typical use for this program will be to convert a gencart output
to gff3 format as:

  cnv_gencartdump2gff.pl -i gencart_dump.txt -o result.gff3

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not require a configuration file or rely on
options set in the user environment.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=over

=item GenCart

The gencart package is available at:


=back

=head1 BUGS AND LIMITATIONS

Known bugs and limitations will be listed here.

=head1 REFERENCE

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

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
