#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_trnascan2gff.pl - Convert trnascan result to GFF3     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/08/2013                                       |
# UPDATED: 01/08/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert primary output from tRNAScan-SE to GFF3 fromat   |
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

# Vars to be set below and used to generate GFF file output
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

# Vars that can be changed to modify GFF output
my $source = "tRNAScan-SE";
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

#-----------------------------+
# FILE HANDLES                |
#-----------------------------+

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

print GFFOUT "##gff-version 3\n";

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

#Sequence                                tRNA    Bounds  tRNA    Anti    Intron Bounds   Cove
#Name                            tRNA #  Begin   End     Type    Codon   Begin   End     Score
#--------                        ------  ----    ------  ----    -----   -----   ----    ------

my $line_num = 0;
while (<INFILE>) {

    # Reset has intron to null
    my $has_intron = 0;

    chomp;
    $line_num++;
    # This skips the header
    next if $line_num < 4;

    print "#IN:".$_."\n" if $verbose;

    my @tr_parts = split( /\s+/, $_);
    $seq_id = $tr_parts[0];
    my $b_begin = $tr_parts[2];
    my $b_end = $tr_parts[3];
    $score = $tr_parts[8];
    $aa = $tr_parts[4];
    my $triplet = $tr_parts[5];
    my $int_begin = $tr_parts[6];
    my $int_end = $tr_parts[7];

    # Set start and end of tRNA strand
    # and set start < end
    if ($b_begin < $b_end) {
	$start = $b_begin;
	$end = $b_end;
	$strand = "+";
    }
    elsif ($b_begin > $b_end) {
	$start = $b_end;
	$end = $b_begin;
	$strand = "-";
    }


    # We have tRNA with introns
    # and mod to make start < end
    my $i_start;
    my $i_end;
    my $i_strand;
    if ($int_begin > 0) {
	$has_intron = 1;
	if ($int_begin < $int_end) {
	    $i_strand = "+";
	    $i_start = $int_begin;
	    $i_end = $int_end;
	}
	elsif ($int_begin > $int_end) {
	    $i_strand = "-";
	    $i_start = $int_end;
	    $i_end = $int_begin;
	}
    }
	

#    my $attribute = "ID=".$seq_id."_".$start."_".end
    my $parent_trna = $seq_id."_".$source."_".$start."_".$end;
    my $parent_name = $type."-".$aa."(".$triplet.")";
    my $attribute = "ID=".$parent_trna.";".
#	"Name=".$type."-".$aa."(".$triplet.")";
	"Name=".$parent_name;
    
    print GFFOUT $seq_id."\t".
	$source."\t".
	$type."\t".
	$start."\t".
	$end."\t".
	$score."\t".
	$strand."\t".
	$phase."\t".
	$attribute.
	"\n";

    if ($has_intron) {
	my $int_attribute = "ID=".$parent_trna."_intron;".
	    "Name=".$parent_name."_intron".";".
	    "Parent=".$parent_trna;

	# Leaving score blank since this is not
	# an intron specific score that the parnet
	# sets.
	print GFFOUT $seq_id."\t".
	    $source."\t".
	    "tRNA_intron\t".
	    $i_start."\t".
	    $i_end."\t".
	    ".\t".
	    $strand."\t".
	    $phase."\t".
	    $int_attribute.
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

cnv_trnascan2gff.pl - Convert tRNAScan-SE primary output to GFF3 format

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    cnv_trnscan2gff.pl -i primary_out.txt -o trnascan_trna.gff3

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

Converts the primary output file from tRNAScan-SE to GFF3 style output. 
This has NOT been extensively tested but it does work for me.

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

STARTED: 01/08/2013

UPDATED: 01/08/2013

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
