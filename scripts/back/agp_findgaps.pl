#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_findgaps.pl - Annotate gaps in dir of fasta files   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/01/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
# Modofied from DAWGPAWS program to find gaps. The AGP      |
# does not appear to have all gaps appropriately            |
# annotated.                                                |
#                                                           |
# LICENSE:                                                  | 
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |
#                                                           |
#-----------------------------------------------------------+
#
# TO DO: Must find a way to speed this up
#

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
#use strict;
use File::Copy;
use Getopt::Long;              # Get options from the command line
use Bio::SeqIO;                # Allows for treatment of seqs as objects
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # To convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 699 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# VARS WITH DEFAULT VALUES
# The following to variables may be leftover cruft 09/29/2008
my $out_ext = ".hard.fasta";  # Outfile extension
my $mask_char = "N";          # Character to mask with

my $ap_path = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";

# BOOLEANS
my $do_apollo_game = 0;        # Convert the gff output to game file format
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $is_gap = 0;                # Current character matche the gap char
my $prev_gap = 0;              # Previous character matches the gap char

# PACKAGE LEVEL SCOPE
my $file_to_mask;              # Path to the file to mask
my $hard_mask_out;             # Path to the hardmasked file output
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc
my $min_gap_len = "100";       # Minimum length to be considered a gap
my $gap_char = "N";            # Character indicating a gap
my $dir_game_out;              # Dir to hold the game xml output
my $dir_gff_out;               # Dir to hold the gff output
my $outfile;
my $infile;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    # Optional strings
		    "gff-ver=s"    => \$gff_ver,
		    "l|len=i"      => \$min_gap_len,
		    "c|gapchar=s"  => \$gap_char,
		    "logfile=s"    => \$logfile,
		    "ap-path=s"    => \$ap_path,
		    # Booleans
		    "game"         => \$do_apollo_game,
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);


my $bac_parent_dir = $outdir;  

my ( $ind_lib , $RepMaskCmd, $MakeGffDbCmd, $MakeGffElCmd );
my ( @RepLibs );
my $ProcNum = 0;

#//////////////////////
my $file_num_max = 2;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

#-----------------------------+
# STANDARDIZE GFF VERSION     |
#-----------------------------+
unless ($gff_ver =~ "GFF3" || 
	$gff_ver =~ "GFF2") {
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
    print "\nbatch_findaps.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">>$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_findgaps.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
}



#-----------------------------+
# CHECK FOR GAPS FOR EACH     |
# SEQUENCE                    |
#-----------------------------+

# OPEN THE SEQUECE OBJECT
my $inseq = Bio::SeqIO->new(-file   => "<$infile",
			    -format => 'fasta' ) 
    || die "ERROR Can not open infile:\n $fasta_file\n";


#-----------------------------+
# OPEN THE OUTPUT FILE        |
#-----------------------------+
if ($outfile) {
    open (GFFOUT, ">$outfile") ||
	die "can not open outfile $outfile\n";
}
else {
    open (GFFOUT, ">&STDOUT") ||
	die "can not write output to STDOUT\n";
}

# Print gff-version header for GFF3 output
if ($gff_ver =~ "GFF3") {
    print GFFOUT "##gff-version 3\n";
}

#-----------------------------+
# CHECK FOR GAPS IN THE OBJ   |
#-----------------------------+
my $time_start_seq = time;
my $attribute;

while (my $seq = $inseq->next_seq()) {
    
    my $seq_len = $seq->length();
    
    print STDERR "FULL SEQ LEN:\t$seq_len\n\n" if $verbose;
    
    my $seq_string = uc($seq->seq());
    # FIND ALL STRINGS THAT ARE RUNS OF N AT LEAST 100 CHARACTERS LONG
    my $gap_count=0;
#    while ($seq_string =~ m/(N{100,})/g ) {
# The following works
#    while ($seq_string =~ m/(N{$min_gap_len,})/g ) {
#    while ($seq_string =~ m/(N{$min_gap_len,})/g ) {
    while ($seq_string =~ m/(($gap_char){$min_gap_len,})/g ) {

	$gap_count++;
	my $gap_len = length($1);
	my $end =  pos($seq_string);
	my $start = $end - $gap_len + 1;
	#print "Found it with end at postion $start \- $end \n";
	#print "\t length: $gap_len\n";
	#print "\tcheck:".$seq->subseq($start,$end)."\n";

	$attribute =  $seq->primary_id."_gap_len_".
	    $gap_len."-".
	    $gap_count;

	print GFFOUT $seq->primary_id()."\t". # Seqname
	    "dawgpaws_findgap\t".             # source
	    "gap\t".                 # feature
	    "$start\t".                       # Start
	    "$end\t".                         # End
	    "$gap_len\t".                     # Score
	    "+\t".                            # Strand
	    ".\t".                            # Frame
	    $attribute."\n";
	
    }


    
} # End of while next_seq


close GFFOUT;

close LOG if $logfile;

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

sub apollo_convert {
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my $InFile = $_[0];        # Input file path
    my $InForm = $_[1];        # Output file format:
                               # game|gff|gb|chadoxml|backup
    my $OutFile = $_[2];       # Output file path
    my $OutForm = $_[3];       # Ouput file foramt
                               # chadoDB|game|chadoxml|genbank|gff|backup
    my $SeqFile = $_[4];       # The path of the sequence file
                               # This is only required for GFF foramt files
                               # When not required this can be passed as na
    my $DbPass = $_[5];        # Database password for logging on to the 
                               # chado database for reading or writing.
    my ( $ApPath, $ApCmd );


# The following path is for the jlb10 machine
#    $ApPath = "/home/jestill/Apps/Apollo/Apollo";
#    $ApPath = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";
#    $ApPath = "/Applications/Apollo/bin/apollo";

    $ApPath = $ap_path;

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ApPath." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}

1;
__END__

=head1 NAME

batch_findgaps.pl - Annotate gaps in a fasta file

=head1 VERSION

This documentation refers to batch_findgaps version $Rev: 699 $

=head1 SYNOPSIS

=head2 Usage

    batch_findgaps.pl -i DirToProcess -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Runs the RepeatMasker program for a set of input
FASTA files against a set of repeat library files &
then converts the repeat masker *.out file into the
GFF format and then to the game XML format for
visualization by the Apollo genome anotation program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

=item --gapchar

The character that is treated as the gap character. B
By default this is N. This option takes a single character 
as its argument.

=item --len

The minimum gap length. This option takes a integer as its option.

=item --game

Use this option to generate a game xml file of the output. This option
requires that you have apollo on your local machine since the program
uses apollo to translate from gff to game xml.

=item --ap-path

Specify the path to your local installation of apollo.

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

=item --test

Run the program without doing the system commands.

=back

=head1 DIAGNOSTICS

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

No configuration files or environmental variables are required to use
this program.

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item Apollo

This program requires the Apollo Genome Annotation Curation tool to
convert the gff output to the game.xml format. This can be obtained at
http://apollo.berkeleybop.org/current/index.html. While Apollo is used to 
convert to game.xml, the batch_findgaps.pl program can be used without
apollo to generate gff foramt files.

=back

=head2 Required Perl Modules

=over 2

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * Incorrect end position of annotated gap

Early versions of this script would assign an incorrect end position the
the 3' end of a gap.

=item * Bug Reporting

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Recognized gap characters

Due to the way that regular expressions are coded in PERL, the characters
that can be used to indicate gaps must be hard coded. The characters that are
currently hard coded for recognition by batch_findgaps are n, N, x, and X. If 
there are additional characters you would like to add as a recognized
gap character, file a Feature Request on the DAWG-PAWS poject page on
Sourceforge ( http://sourceforge.net/tracker/?group_id=204962&atid=991722 ).

=item * Limited file extensions are supported

BLAST output file must currently end with blo, bln, or blx. For example
a BLASTx output may be named BlastOut.blx while a BLASTN output
may be names BlastOut.bln. FASTA files must end with a fasta or fa extension.
For examples must have names like my_seq.fasta or my_seq.fa.

=back

=head1 SEE ALSO

The batch_findgaps.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/01/2007

UPDATED: 09/29/2008

VERSION: $Rev: 699 $

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 08/01/2007
# - Program started
# - Base program and POD docs written
#
# 12/06/2007
# - Moved POD documentation to the end of the file
# - Added info to POD documentation
# 
# 05/20/2008
# - Added some optional output for the verbose mode
# - Removed commented out code no longer needed
# - The current implementation is too slow for anything
#   larger then a BAC sized contig
#     - Using seq->subseq
#     - For a contig of size 7,822,695 this would 
#       take about 1,431 hours to run
#      --may be quicker to just load seq string to an array and pop
#
# 09/29/2008
# - Made the conversion to game xml format with apollo a boolean
# - Updated help information
# - Added option to specify the Apollo path for convertion to
#   game xml. This is a sloppy implementation and uses a global variable 
#   in the subfunction.
