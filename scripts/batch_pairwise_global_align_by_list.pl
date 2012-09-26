#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_pairwise_global_align.pl
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/07/2012                                       |
# UPDATED: 09/24/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Do needle alignment and report values that exceed        |
#  threshold. This takes a tab delimited input file         |
#  to determine the pairs of sequences to do the            |
# aliggments for. THis will allow for a step                |
# where global alignments are only sought for pairs         |
# that have a signifiicant BLAST alignment.                 |
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

#use Bio::SimpleAlign;          # Use bioperl simple alignment to parse 
use Bio::AlignIO; 
use Bio::SeqIO;                # use to get the sequence length


print STDERR "Working\n";

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/ || "pre-release";
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;           # This is the file listing the files to search
my $outdir;
my $matrix_out;
my $indir;            # This is the source of the fasta files to process

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $do_self = 0;
my $do_keep_temp = 0;             # Keep the temp file

my $a_file_path;
my $b_file_path;

# The threoshold value to report in the matrix
my $threshold = 80;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"     => \$indir,
		    "f|infile=s"    => \$infile,
                    "o|outdir=s"    => \$outdir,
		    # outfile for matrix output
		    "s|self"        => \$do_self,
		    "m|matrix=s"    => \$matrix_out,
		    # ADDITIONAL OPTIONS
		    "k|keep-temp"   => \$do_keep_temp,
		    "t|threshold=s" => \$threshold,
		    "gff-ver=s"     => \$gff_ver,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "test"          => \$do_test,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

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

if ($matrix_out) {
    open (MATOUT, ">$matrix_out") ||
	die "Can not open matrix output file\n";
}
else {
    open (MATOUT, ">&STDOUT") ||
	die "Can not open STDOUT for writing\n";
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

# Get all of the fasta files in the input dir

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
#unless ($indir =~ /\/$/ ) {
#    $indir = $indir."/";
#}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

# Open input file
if ($infile) {
    open (INFILE, "<$infile") ||
	die "Can not open intput file\n";
    print STDERR "Processing file $infile\n";
}
else {
    die "An input file is required\n";
}

#-----------------------------+
# CREATE OUTPUT DIR IF NEEDED
#-----------------------------+
unless (-e $outdir) {
    mkdir $outdir ||
	die "Could not create dir:\n$outdir\n"
}


#-----------------------------+
# DO AN ALIGNMENT             |
# FOR EACH SEQUENCE IN THE    |
# INPUT LIST                  |
#-----------------------------+
my $file_num = 0;
my $proc_num = 0;
my $line_num = 0;
while (<INFILE>) {
    chomp;
    $line_num++;

    # Allow for comment lines preceeding by #
    next if m/^\#/; 

    print STDERR "$line_num\n";

    my @in_parts = split(/\t/, $_);
    my $num_parts = @in_parts;
    # Only if we have qry and qry
    # Allowing for 3 cols lets me use the sim file as input
    if ($num_parts == 2 ||
	$num_parts == 3) {
	
	my $a_name_root = $in_parts[0];
	my $b_name_root = $in_parts[1];
#	my $a_path = $indir.$a_name_root.".fasta";
#	my $b_path = $indir.$b_name_root.".fasta";

	my $a_path = $indir.$in_parts[0].".fasta";
	my $b_path = $indir.$in_parts[1].".fasta";
#	my $b_path = $indir.$in_parts[1].".fasta";
     
	print STDERR "A in : ".$in_parts[0]."\t";
	print STDERR "B in : ".$in_parts[1]."\n";

	print STDERR "A path: ".$a_path."\n";
	print STDERR "B path: ".$b_path."\n";

	# Need to get the length of the query sequence to
	# be able to do the comparision.
	my $a_seq_in  = Bio::SeqIO->new(-file => $a_path , 
					'-format' => 'fasta');
	my $a_seq_obj = $a_seq_in->next_seq();
	my $a_seq_length = $a_seq_obj->length;

	my $ab_out = $outdir.$a_name_root."_".$b_name_root.".needle";
	
	my $cmd = "needle".
	    " -asequence ".$a_path.
	    " -bsequence ".$b_path.
	    " -outfile ".$ab_out.
	    " -gapextend 0.5".
	    " -gapopen 10 ".
	    " 2>/dev/null";

	print STDERR "$cmd\n" 
	    if $verbose;

	system ($cmd) unless $do_test;
		
	open (RESULT , "<$ab_out" ) ||
	    die "Can not open the result file\n";

	# Just get what we need from the result file 
	# which is the identity
	my $num;
	my $denom;
	my $percent_id;
	while (<RESULT>) {
	    chomp;
	    if (m/^# Identity:\s*(.*)\/(.*) \((.*)\)/) {
		$num = $1;
		$denom = $2;
		$percent_id = $3;
		last;
	    }
	}

	my $a_percent_align = ($num/$a_seq_length) * 100;
	$a_percent_align = sprintf("%.2f", $a_percent_align);

	# Look at hte reuslt

	if ($verbose) {
	    print STDERR "\n";
	    print STDERR "A-seq length: ".$a_seq_length."\n";
	    print STDERR "aln_num: ".$num."\n";
	    print STDERR "aln_denom: ".$denom."\n";
	    print STDERR "aln_percent: ".$percent_id."\n";
	    print STDERR "a_percent_algn: ".$a_percent_align."\n";
	}
	
	if ($a_percent_align >= $threshold) {
	    print STDERR "MATCH\n" if $verbose;
	    print MATOUT $a_name_root."\t".$b_name_root."\t".
		$a_percent_align."\n";
	}
	
	# Remove the temp files unless we want to keep them
	unless ($do_keep_temp) {
	    unlink ($ab_out);
	}

    }

}

close (INFILE);

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

batch_pairwise_global_align.pl - Needleman–Wunsch in batch mode

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    batch_pairwise_global_align.pl -i indir -o outdir -m matrix.txt \
                                   -t 80

=head2 Required Arguments

    --indir,i          # Path to the output dir
    --outdir,-o        # Path to the input dir
    --threshold,-t     # Minimum threshold

=head1 DESCRIPTION

Run pairwise global alignment for all possible pairwise comparisons
given a directory of fasta files where each file is a fasta file
containing a single squence. This will produce a matrix of all pairwise
matches that meet a minimum threshold of identity. For example running 
this program with the deafult threshold of 80 will report all pairwise
sets that meet the 80/80 threshold as set by Wicker et al to 
define transposable element families. This reports the matrix 
as a similarity matrix suitable for use with the repminer program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path to the input dir containing the fasta files to align.

=item -o,--outfile

Path to the output dir where the results will be stored.

=back

=head1 OPTIONS

=over 2

=item -s,--self

Include self comparisons. By default, self comparisoins are not included.

=item -k,--keep-temp

This will keep the results of the needleman alignment.

=item -t,--threshold

The minimumn threshold value to report. This must be a number between
0 and 100. Be default this value is 80.

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

To run the program with default settings

 batch_pairwise_global_align.pl -i input_dir -o output_dir -m matrix_out.txt

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

=head2 Software

=over

=item EMBOSS::needle

This program requires the needle global alignment program from EMBOSS.
http://emboss.sourceforge.net/

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

Please refer to the DAWGPAWS manuscript in Plant Methods when describing
your use of this program:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/07/2012

UPDATED: 09/07/2012

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
