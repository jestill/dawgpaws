#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_aragorn2gff.pl - Convert aragorn trna output to gff   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/08/2013                                       |
# UPDATED: 01/08/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert a aragorn tab delimited output from the          |
#  aragorn tRNA prediction program to a usable GFF file.    |
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

# Vars that can be changed to modify GFF output
my $source = "Aragorn";
my $type = "tRNA";


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "source=s"    => \$source,
		    "type=s"      => \$type,
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
	die "Can not open file for input at $infile\n";
}

if ($outfile) {
    open (GFFOUT, ">$outfile") ||
	die "Can not open output file for output at $outfile";
}
else {
    open (GFFOUT, ">&STDOUT" ) ||
	die "Can not open STDOUT for output. Specify outfile with -o option\n"
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

my $seq_id;
my $start;
my $end;
my $score = ".";
my $strand;
my $phase = ".";
my $attribute;
my $gff_str;
my $num_result = 0;
my $aa;

while (<INFILE>) {

    chomp;
#    print STDERR $_."\n";

    if (m/^\>(.*)/) {
#	print STDERR $_."\n";
	$seq_id = $1;
#	print STDERR "\t".$seq_id."\n";
	
    }
    elsif (m/^(.*)\t(.*)\t(.*)/) {
	$num_result++;

	my ($num,$trna_type,$loc) = split(/\s+/ , $1);
	my $len = $2;
	my $triplet = $3;
	# Make triplet code uppercase
	$triplet = uc($triplet);


	if ($verbose) {
	    print STDERR "\t\tSeq: ".$seq_id."\n";
	    print STDERR "\t\tNum: ".$num."\n";
	    print STDERR "\t\tType: ".$trna_type."\n";
	    print STDERR "\t\tLoc: ".$loc."\n";
	    print STDERR "\t\t$len\n";
	    print STDERR "\t\t$triplet\n";
	}

	# Type will be used for Name=
	# Will need to get start and end below
#	if ($loc =~ m/c.*/ ) {
	if ($loc =~ m/c\[(.*)\,(.*)\]/ ) {
#	    # I think this is complement strand
	    $start = $1;
	    $end = $2;
	    $strand = "-";
	}
	elsif ($loc =~ m/\[(.*)\,(.*)\]/ ) {
	    $start = $1;
	    $end = $2;
	    $strand = "+";
	}
	else {
	    $strand = "+";
	    $start = "ERR";
	    $end = "ERR";
	}


	if ($trna_type =~ m/(.*)\-(.*)/) {
	    $type = $1;
	    $aa = $2;
	}


	my $attribute = "ID=".$seq_id."_".$source."_".$start."_".$end.";".
#	    "Name=".$trna_type;
	    "Name=".$type."-".$aa.$triplet;

	# TEST GFF 3 STRING
	print GFFOUT $seq_id."\t".
	    $source."\t".
	    $type."\t".
	    $start."\t".
	    $end."\t".
	    $score."\t".
	    $strand."\t".
	    $phase."\t".
	    $attribute."\t".
	    "\n";
	
	
    }

  

}


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

cnv_aragorn2gff.pl

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    cnv_aragorn2gff.pl -i aragorn_result.txt -o ar_out.gff

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

Converts a basic aragorn tab delimited text file to a GFF3 format
result. This currently does not support GFF2 output. This has been tested
with aragorn result that was run as:

 ./aragorn -rp -w -seq -l -t ambo.fasta -o ara_test_out.txt

Currently will tag all results as tRNA.

using Aragorn v.1.2.34

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

This is a typcial use case.

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

= head2 Software

This program required the Aragorn tRNA prediction program:
http://mbio-serv2.mbioekol.lu.se/ARAGORN/

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

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

STARTED: 01/08/2013

UPDATED: 01/08/2013

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
