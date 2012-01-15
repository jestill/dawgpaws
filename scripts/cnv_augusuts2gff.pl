#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_augustus2gff.pl - Convert augustus output to 
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/15/2012                                       |
# UPDATED: 01/15/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  While Augustus does produce output in a format that can  |
#  considered a variant of GFF3, this parses the Augustus   |
#  GFF to a more standard GFF3 format that is used by       |
#  DAWGPAWS. The goal is to produce a Name and ID value     |
#  that is globally unique within the genome annotation     |
#  set.                                                     |
#                                                           |
# USAGE:                                                    |
#  cnv_augustus2gff.pl -i aug_output -o gff_output          |
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

my $prog_source = "Augustus";
my $prog_suffix = 0;
my $delim = "_";
my $sequence_id = 0;


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

print STDERR "GFF VERSION is $gff_ver\n" if $verbose;

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
    print "\ncnv_augustus2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
# Take the input file and convert to GFF3
augustus2gff ($prog_source, $infile, $outfile, $sequence_id, $prog_suffix, 
	      $delim, 0, 0);


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


sub augustus2gff {
# Technically Augustus returns data in GFF format, this modifies the 
# Augustus native output to the GFF format that is standard for 
# DAWGPAWS

    
    my ($source, $aug_in, $gffout, $seq_id, $src_suffix, $d,
	$do_append, $do_apollo_format) = @_;

    my $i = -1;
    my $j = -1;
    my $model_num = 0;

    # The following delimiter (:) works for Apollo, but throws problems
    # with GBROWSE so switching to a variable delim
    if ($src_suffix) {
	$source = $source.$d.$src_suffix;
    }

    my $attribute;

    #-----------------------------+
    # OPEN INPUT FILE HANDLE      |
    #-----------------------------+
    if ($aug_in) {
	open (INFILE, "<$aug_in") ||
	    die "ERROR: Can not open Augustus result file\n $aug_in\n";
	
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }
    
    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gffout) {
	if ($do_append) {
	    open (GFFOUT, ">>$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
	else {
	    open (GFFOUT,">$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }


    #-----------------------------+
    # PROCESS AUGUSTUS OUTPUT
    #-----------------------------+
    my $prev_gene_name = "NULL";
    my $exon_count = 0;            # Counter for exons within gene
    my $intron_count = 0;          # Counter for introns within gene
    my $cds_count = 0;
    while (<INFILE>) {
	chomp;
	# Skip comment lines
	next if m/^\#/;

	print STDERR $_."\n"
	    if $verbose;

	my @gff_parts = split;
	my $num_gff_parts = @gff_parts;

	# Reset exon and intron counters, these counters
	# allow for unique ids to be assigned to introns
	# exons 
	if ($gff_parts[2] =~ "gene") {
	    $exon_count=0;
	    $intron_count=0;
	    $cds_count=0;
	}

	#-----------------------------+
	# MAKE START < END            |
	#-----------------------------+
	my $start;
	my $end;
	if ($gff_parts[3] < $gff_parts[4]) {
	    $start = $gff_parts[3];
	    $end = $gff_parts[4];
	} else {
	    $end = $gff_parts[4];
	    $start = $gff_parts[3];
	}

	my $type = $gff_parts[2];
	my $attribute = $gff_parts[8];

	# WE MAY NEED TO ORDER PARTS OF A SINGLE GENE
	# SEQUENTIALY BY OCCURENCE ON THE CONTIG

	# It may also be necessary to include contig name in result
	# name to that results for multiple contigs can have
	# unique IDs.
	# Load results to the $aug_results hash
#	$i++;
#	$aug_results[$i][0] = $gff_parts[0];
#	$aug_results[$i][1] = $program.":".$src_suffix;
#	$aug_results[$i][2] = $gff_parts[2];
#	$aug_results[$i][3] = $start;
#	$aug_results[$i][4] = $end;
#	$aug_results[$i][5] = $gff_parts[5];
#	$aug_results[$i][6] = $gff_parts[6];
#	$aug_results[$i][7] = $gff_parts[7];
#	$aug_results[$i][8] = $gff_parts[8];


	#/////////////////
	# Kludge alert
	#/////////////////
	# For GFF2 will only report the CDS
	if ($gff_ver =~ "GFF2") {
	    next unless $type =~ "CDS";

	    # Replace longer GFF3 type name
	    if ($attribute =~ m/Parent\=(.*)/) {
		$attribute = $1;
	    }
	    
	    # Apollo can't render CDS features so 
	    # converting type to exon for Apollo
	    if ($do_apollo_format) {
		$type = "exon";
	    }
	}
	else {
	    # Add seuqence name to the ID string to make globally
	    # unique within the context of the assembly
	    $attribute =~ s/ID=/ID=$gff_parts[0]./g;
	    $attribute =~ s/Parent=/Parent=$gff_parts[0]./g;

	    # Add uniquie ID to cds
#	    if ( $gff_parts) {
	    if ( $gff_parts[2] =~ "CDS") {
		$cds_count++;
		$attribute =~ s/cds/cds.$cds_count/g;
#		print STDERR "CDS\n";
	    }

	    # Add an id string to exons and introns
	    if ($attribute =~ m/^Parent=(.*)/ ) {
		my $id_root = $1;
		if ( $gff_parts[2] =~ "exon") {
		    $exon_count++;
		    $attribute = "ID=".$id_root.".ex.".$exon_count.";"
			.$attribute
		}
		elsif ( $gff_parts[2] =~ "intron") {
		    $intron_count++;
		    $attribute = "ID=".$id_root.".in.".$intron_count.";"
			.$attribute;
		}
		elsif ( $gff_parts[2] =~ "CDS") {
		    $cds_count++;
		    $attribute = "ID=".$id_root.$cds_count.";"
			.$attribute;
		}
		elsif ( $gff_parts[2] =~ "transcription_start_site") {
		    $attribute = "ID=".$id_root.".tss;"
			.$attribute;
		    $type = "TSS";
		}
		elsif ( $gff_parts[2] =~ "start_codon") {
		    $attribute = "ID=".$id_root.".start;"
			.$attribute;
		}
		elsif ( $gff_parts[2] =~ "stop_codon") {
		    $attribute = "ID=".$id_root.".stop;"
			.$attribute;
		}
	    }


	    # Get name from ID and append to end of the attribute
	    # string for genes
	    if ( $gff_parts[2] =~ "gene" ) {
		if ($attribute =~ m/ID=(.*);/) {
		    $attribute = $attribute.";Name=".$source.".".$1;
		}
		elsif ($attribute =~ m/ID=(.*)/) {
		    $attribute = $attribute.";Name=".$source.".".$1;
		}
	    }

	}

	# PRINT GFF OUTPUT
	print GFFOUT $gff_parts[0]."\t".
	    $source."\t".
	    $type."\t".
	    $start."\t".
	    $end."\t".
	    $gff_parts[5]."\t".
	    $gff_parts[6]."\t".
	    $gff_parts[7]."\t".
	    $attribute."\n";

    } # End of while INFILE

    close INFILE;
    close GFFOUT;

}


1;
__END__

=head1 NAME

cnv_augustus2gff.pl - Converts Augustus output to GFF3 standard

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
