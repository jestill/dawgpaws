#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fgenesh2gff.pl - Convert fgenesh output to gff        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/31/2009                                       |
# UPDATED: 01/31/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert output from fgenesh to the gff format.           |
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
use Bio::Tools::Fgenesh;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $seqname;
my $param;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $append = 0;
my $test = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "p|param=s"   => \$param,
		    "n|name=s"    => \$seqname,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "append"      => \$append,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

# MAY NEED TO ITERATE ACROSS THE FILE AND GET RID OF COPYWRITE STATEMENT

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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# DO THE CONVERSION           |
#-----------------------------+
my $prog = "fgenesh";
fgenesh2gff ($prog, $infile, $outfile, $seqname, $param, $append);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub fgenesh2gff {
    # fgnesh_in  - path to the fgenesh program
    # gff_out    - path to the gff output file
    # seq_id     - id of the source sequence
    # src_suffix - parameter id for fgenesh run

    my ($source, $fgenesh_in, $gff_out, $seq_id, $src_suffix, $do_append ) = @_;

    #-----------------------------+
    # OPEN THE FGENESH INFILE     |
    #-----------------------------+
    my $fgenesh_result;
    if ($fgenesh_in) {
	$fgenesh_result = Bio::Tools::Fgenesh->new(-file => $fgenesh_in);
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	$fgenesh_result = Bio::Tools::Fgenesh->new( -fh  => \*STDIN );
    }

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
     # Default to STDOUT if no arguemtn given
    if ($gff_out) {
	if ($do_append) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
	else {
	    open (GFFOUT,">$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    #-----------------------------+
    # SET PROGRAM SOURCE          |
    #-----------------------------+
    unless ($source) {
	$source = "fgenesh";
    }
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    my $gene_num = 0;
    while (my $gene = $fgenesh_result->next_prediction()) {

	$gene_num++;
	#-----------------------------+
	# SET SEQUENCE ID             |
	#-----------------------------+
	unless ($seq_id) {
	    if ($gene->seq_id()) {
		$seq_id = $gene->seq_id();
	    }
	    else {
		$seq_id = "seq";
	    }
	}

	# $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
	# off Bio::SeqFeature::Gene::Transcript.
	#
	# $gene->exons() returns an array of 
	# Bio::Tools::Prediction::Exon objects
	# all exons:
	my @exon_arr = $gene->exons();
	

	foreach my $ind_exon (@exon_arr) {
	    #print STDERR $ind_exon;
#	    if ($ind_exon->is_coding()) {
#		print STDERR "coding\t";
#	    }
#	    else {
#	    }

	    #-----------------------------+
	    # FORMAT STRAND               |
	    #-----------------------------+
	    my $strand = $ind_exon->strand()."\t";
	    if ($strand =~ "-1") {
		$strand = "-";
	    }
	    else {
		$strand = "+";
	    }

	    #-----------------------------+
	    # GET START AND END           |
	    #-----------------------------+
	    my $start = $ind_exon->start();
	    my $end = $ind_exon->end();

	    if ($start > $end) {
		$end =  $ind_exon->start();
		$start = $ind_exon->end();
	    }

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    #-----------------------------+
	    # The following prints one line at a time
	    #print GFFOUT $seq_id."\t";     # Seqname
	    #print GFFOUT $source."\t";     # Source
	    #print GFFOUT "exon\t";         #feature
	    #print GFFOUT $start."\t";      # start
	    #print GFFOUT $end."\t";        # end
	    #print GFFOUT $ind_exon->score()."\t";  # score
	    #print GFFOUT $strand."\t";       # strand
	    #print GFFOUT ".\t";              # frame
	    #print GFFOUT "gene_".$gene_num;  # attribute
	    #print GFFOUT "\n";

	    print GFFOUT $seq_id."\t".     # Seqname
		$source."\t".     # Source
		"exon\t".         #feature
		$start."\t".      # start
		$end."\t".        # end
		$ind_exon->score()."\t".  # score
		$strand."\t".       # strand
		".\t".              # frame
		"gene_".$gene_num.  # attribute
		"\n";


	    # The following does work
	    #print GFFOUT $ind_exon->primary_tag()."\t";


	    #///////////////////////////////
	    # The following do not work
	    #///////////////////////////////
	    #print GFFOUT $gene->cds()."\n";
	    #print GFFOUT $ind_exon->significance()."\t";
	    #print $ind_exon->predicted_cds();
	    #print GFFOUT $ind_exon->coding_signal_score()."\t";
	    #print GFFOUT $ind_exon->seq_id()."\t";
	    #print $ind_exon->significance()."\n";
	    # Get the CDS of the sequence
	    #print STDERR $ind_exon->cds()."\n";

	}
	
#       # initial exons only
#       @init_exons = $gene->exons('Initial');
#       # internal exons only
#       @intrl_exons = $gene->exons('Internal');
#       # terminal exons only
#       @term_exons = $gene->exons('Terminal');
#       # singleton exons: 
#       ($single_exon) = $gene->exons();


   }

    # CLOSE FGENESH
    $fgenesh_result->close();
    

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InFile -o OutFile

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

=head2 Required Software

=head2 Required Perl Modules

=over 2

= * Bio::Tools::Fgenesh

This program requires the perl module Bio::Tools::Fgenesh. This module is
part of the bioperl package

=back

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

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
