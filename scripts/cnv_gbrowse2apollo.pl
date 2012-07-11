#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gbrowse2apollo.pl - Convert GBrowse GFF to Apollo GFF |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/29/2012                                       |
# UPDATED: 07/01/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert GFF3 exported from GBrowse to a format that is   |
#  compatible with the way the Apollo treats GFF files.     |
#                                                           |
# USAGE:                                                    |
#  cnv_gbrowse2apollo.pl -i infile.gff3 -o outfile.gff3     |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# Curretly this is limited to taking GFF3 GBrowse data as input and producing
# GFF3 Apollo data as output.

package DAWGPAWS;

# To do
# Strip out URL encoding of Target
# %3D should be =
# %20 should be space

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
my $gff_ver = uc($ENV{DP_GFF}) || "GFF3";

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

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
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

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

if ($infile) {
    open(INFILE, "<$infile") ||
	die "Can not open infile $infile";
}
else {
    print STDERR "Expecting input from STDIN\n";
    open (INFILE, "<&STDIN") ||
	die "Can not open STDIN for input";
}

if ($outfile) {
    open(OUTFILE, ">$outfile") ||
	die "Can not open output file $outfile";
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not open STDOUT for output";
}

# Write GFF3 header to output file
print OUTFILE "##gff-version 3\n";

# Pragma star
my $seq_region_id;
my $seq_region_start;
my $seq_region_end;
my $in_seq_data = 0;
my $in_gff_data = 0;
while (<INFILE>) {

    chomp;

    # Get pragmas and determine if we are in sequence data region
    if (m /^#/ ) {
	
	$in_gff_data = 1;
	
	if (m /^##sequence-region\s(.*)\s(.*)\s(.*)/) {
	    print STDERR "In sequence region pragma\n"
		if $verbose;
	    $seq_region_id = $1;
	    $seq_region_start = int($2);
	    $seq_region_end = int($3);
	    print STDERR "\t".$seq_region_id."\t".
		$seq_region_start."\t".
		$seq_region_end."\n"
		if $verbose;
	    print OUTFILE $_."\n";

	}
	# The following appears to be the nature of the 
	# gff3 header in some cases. Perhaps when the
	# sequence file is not included in the output
	elsif (m /^##sequence-region\s(.*)\:(.*)\.\.(.*)/) {
	    $seq_region_id = $1;
	    $seq_region_start = int($2);
	    $seq_region_end = int($3);
	    print STDERR "\t".$seq_region_id."\t".
		$seq_region_start."\t".
		$seq_region_end."\n"
		if $verbose;
	    print OUTFILE $_."\n";
	}
	# All other comments and pragmas are ignored
	next;
    }
    elsif (m /^>/) {
	$in_seq_data = 1;
	print STDERR "In seq data\n" 
	    if $verbose;
    }


    # Data in the sequence region are printed out without modification
    # GFF tab delim lines may need modification
    if ($in_seq_data) {
	print OUTFILE $_."\n";
    }
    elsif ($in_gff_data) {

	
	my ($seq_id, $source, $type,
	    $start, $end, $score,
	    $strand, $phase, $attribute) = split ('\t', $_);


	#-----------------------------+
	# TO DO: CLEAN ATTRIBUTE HERE OF SPECIAL CHARACTERS
	# USING 

	# Skip stop_codons and start_codons since Apollo does not
	# support these features in GFF3
	if ( $type =~ "stop_codon" ||
	     $type =~ "start_codon" ) {
	    next;
	}

	# skip coverage features
	if ($source =~ "coverage") {
	    next;
	}

	if ($type =~ "coverage") {
	    next;
	}

	# skip sam/bam features
	if ($source =~ "sam/bam") {
	    next;
	}

	# check start stop relationship
	$start = int($start);
	$end = int($end);
	
	# Currently lines falling outside of bounds will be ignored
	# This will probably cause problems with parent child features ...
	# It looks like Apollo has no trouble handeling features that
	# are outside of the scoe of the region.
	if ($start < $seq_region_start ) {
	    print STDERR "Feature start less than region start:\n".
		"$_\n"
		if $verbose;
	    print STDERR $start."<".$seq_region_start."\n"
		if $verbose;
	    next;
	}
	elsif ($end > $seq_region_end) {
	    print STDERR "Feature end greater than region end:\n".
		"$_\n"
		if $verbose;
	    print STDERR $end.">".$seq_region_end."\n"
		if $verbose;
	}

	# For match features
	# print out as they are in the incoming GFF file
	if ($type =~ "match") {
	    print OUTFILE $_."\n";
	    next;
	}
	# Change tandem_repeat feature to match feature
	elsif ($type =~ "tandem_repeat") {
	    print OUTFILE $seq_id."\t".
		$source."\t".
		"match\t".
		$start."\t".
		$end."\t".
		$score."\t".
		$strand."\t".
		$phase."\t".
		$attribute.
		"\n";
	    next;
	}
	# we may be able to leave gap unmodified but this has not
	# been tested to making htis a match feature for now
	elsif ($type =~ "gap") {
	    print OUTFILE $seq_id."\t".
		$source."\t".
		"match\t".
		$start."\t".
		$end."\t".
		$score."\t".
		$strand."\t".
		$phase."\t".
		$attribute.
		"\n";
	    next;
	}


	print OUTFILE $_."\n";
	
    }

#    print STDERR $_."\n";

}

unless ($in_seq_data) {
    print "\a";
    print STDERR "\n\n";
    print STDERR "//////////////////////////////////////////////////////\n";
    print STDERR "//////////////////////////////////////////////////////\n";
    print STDERR "WARING: It appears the sequence file is not included\n".
	" in this GFF file\n";
    print STDERR "//////////////////////////////////////////////////////\n";
    print STDERR "//////////////////////////////////////////////////////\n";
}

close(INFILE);
close(OUTFILE);


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

    Name.pl -i InFile -o OutFile

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

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

Other modules or software that the program is dependent on.

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

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
