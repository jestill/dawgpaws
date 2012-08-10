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
use Bio::SeqIO;               # Read and write seq files
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

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# Global count variables
my $tot_G = 0;
my $tot_C = 0;
my $tot_A = 0;
my $tot_T = 0 ;
my $tot_g = 0;
my $tot_c = 0;
my $tot_a = 0;
my $tot_t = 0;
my $tot_N = 0;
my $tot_n = 0;
my $tot_len = 0;

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
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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

if ($outfile) {
    open (OUTFILE, ">$outfile") ||
	die "Can not open output file\n";
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not STOUT for output\n";
}

my $format = "fasta";
my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
                             '-format' => $format);

print OUTFILE "seq_id\tlen\tratio_masked_len\tratio_masked_bases\t".
    "ratio_n\t".
    "ratio_gc\t\tG\tC\tA\tT\tg\tc\ta\tt\tN\n\n";
while (my $inseq = $seq_in->next_seq) {
    my $id = $inseq->primary_id();

    my $len = $inseq->length;
    $tot_len = $tot_len + $len;

    # Assume uppercase unmasked
    my $G = $inseq->seq() =~ tr/G//; 
    $tot_G = $tot_G + $G;
    my $C = $inseq->seq() =~ tr/C//; 
    $tot_C = $tot_C + $C;
    my $A = $inseq->seq() =~ tr/A//; 
    $tot_A = $tot_A + $A;
    my $T = $inseq->seq() =~ tr/T//; 
    $tot_T = $tot_T + $T;

    # Lower case masked
    my $g = $inseq->seq() =~ tr/g//; 
    $tot_g = $tot_g + $g;
    my $c = $inseq->seq() =~ tr/c//; 
    $tot_c = $tot_c + $c;
    my $a = $inseq->seq() =~ tr/a//; 
    $tot_a = $tot_a + $a;
    my $t = $inseq->seq() =~ tr/t//; 
    $tot_t = $tot_t + $t;

    # N
    my $N = $inseq->seq() =~ tr/N//; 
    $tot_N = $tot_N + $N;
    my $n = $inseq->seq() =~ tr/n//; 
    $tot_n = $tot_n + $n;

    my $percent_mask = sprintf ( "%.4f", ( ($g + $c + $a + $t)/$len) );
    my $percent_mask_real = sprintf ( "%.4f", ( ($g + $c + $a + $t)/( $G + $C + $A + $T + $g + $c + $a + $t )) );
    my $percent_n = sprintf ( "%.4f", (($N + $n)/$len) );
    my $percent_gc = sprintf ( "%.4f", (($G + $g + $C + $c)/$len) );
    my $percent_gc_real = sprintf ( "%.4f", (($G + $g + $C + $c)/ ($G + $C + $A + $T + $g + $c + $a + $t)  ) );
    
    # print output 
    print OUTFILE $id."\t".$len."\t".
	$percent_mask."\t".$percent_mask_real."\t".
	$percent_n."\t".
	$percent_gc."\t".$percent_gc_real."\t".
	$G."\t".$C."\t".$A."\t".$T."\t".
	$g."\t".$c."\t".$a."\t".$t."\t".
	$N."\t".$n.
	"\n";

}

print OUTFILE "-----------------------------------\n";

my $tot_percent_mask = sprintf ( "%.4f", ( ($tot_g + $tot_c + $tot_a + $tot_t)/$tot_len) );
my $tot_percent_mask_real = sprintf ( "%.4f", ( ($tot_g + $tot_c + $tot_a + $tot_t) / ($tot_G + $tot_C + $tot_A + $tot_T + $tot_g + $tot_c + $tot_a + $tot_t ) ) );
my $tot_percent_n = sprintf ( "%.4f", (($tot_N + $tot_n)/$tot_len) );
my $tot_percent_gc = sprintf ( "%.4f", (($tot_G + $tot_g + $tot_C + $tot_c)/$tot_len) );
my $tot_percent_gc_real = sprintf ( "%.4f", ( ($tot_g + $tot_c + $tot_G + $tot_C)/ ($tot_G + $tot_C + $tot_A + $tot_T + $tot_g + $tot_c + $tot_a + $tot_t ) ) );

print OUTFILE "Total\t".$tot_len."\t".
    $tot_percent_mask."\t".$tot_percent_mask_real."\t".
    $tot_percent_n."\t".
    $tot_percent_gc."\t".$tot_percent_gc_real."\t".
    $tot_G."\t".$tot_C."\t".$tot_A."\t".$tot_T."\t".
    $tot_g."\t".$tot_c."\t".$tot_a."\t".$tot_t."\t".
    $tot_N."\t".$tot_n.
    "\n";

# Global output

#if ($infile) {
#    open (INFILE, "<$infile") ||
#	die "Can not open inputfile:$infile";
#}
#else {
#    open (INFILE, "<&STDIN") ||
#	die "can not open stdin for input";
#}


#while (<INFILE>) {
#    chomp;
#}

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
