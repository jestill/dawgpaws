#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_findgaps.pl - Annotate gaps in a fasta file         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/01/2007                                       |
# UPDATED: 12/06/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of fasta files, this will find gaps    |
#  in the fasta file and report these as gaps in a gff file |
#  as well as produce a game.xml file for apollo.           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |
#                                                           |
#-----------------------------------------------------------+
#
# TO DO: Appears to have error in delineating stop of the gap
#
#-----------------------------+
# INCLUDES                    |
#-----------------------------+
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
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# VARS WITH DEFAULT VALUES
my $out_ext = ".hard.fasta";  # Outfile extension
my $mask_char = "N";          # Character to mask with

# BOOLEANS
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
my $gap_char = "n";            # Character indicating a gap
my $dir_game_out;              # Dir to hold the game xml output
my $dir_gff_out;               # Dir to hold the gff output

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "len"          => \$min_gap_len,
		    "gapchar"      => \$gap_char,
		    "logfile=s"    => \$logfile,
		    # Booleans
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
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}


if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\nbatch_mask.pl:\n".
	"Version: $ver\n\n";
    exit;
}

# Show full help when required options
# are not present
if ( (!$indir) || (!$outdir) ) {
    print_help("full");
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
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# CHECK FOR GAPS FOR EACH     |
# SEQUENCE                    |
#-----------------------------+
for my $ind_file (@fasta_files)
{
    
    $file_num++;
    
    # Temp exit for debug
    #if ($file_num == $file_num_max) {exit;}


    print "==============================\n" if $verbose;
    print "  $file_num of $num_files\n" if $verbose;
    print "==============================\n" if $verbose;

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FILE TO MASK                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }

    #-----------------------------+
    # CREATE PARENT DIR IF IT     |
    # DOES NOT ALREADY EXIST      |
    #-----------------------------+
    my $dir_parent = $outdir.$name_root."/";
    unless (-e $dir_parent) {
	print "creating dir: $dir_parent\n";
	mkdir $dir_parent ||
	    die "Could not creat the output dir:\n$dir_parent\n";
    }

    #-----------------------------+
    # MAKE SURE THE GFF AND GAME  |
    # DIRS ARE PRESENT            |
    #-----------------------------+
    $dir_game_out = $outdir.$name_root."/game/";
    unless (-e $dir_game_out) {
	print "Creating dir:\n$dir_game_out\n" if $verbose;
	mkdir $dir_game_out ||
	    die "Could not create dir:\n$dir_game_out\n";
    }
    
    $dir_gff_out = $outdir.$name_root."/gff/";
    unless (-e $dir_gff_out) {
	print "Creating dir:\n$dir_gff_out\n" if $verbose;
	mkdir $dir_gff_out ||
	    die "Could not create dir:\n$dir_gff_out\n";
    }


    #-----------------------------+
    # INPUT AND OUTPUT FILE PATHS |
    #-----------------------------+
    my $fasta_file = $indir.$ind_file;
    my $gff_out = $dir_gff_out.$name_root."_gaps.gff";
    my $game_out = $dir_game_out.$name_root."_gaps.game.xml";
    
    #-----------------------------+
    # LOAD FASTA FILE TO SEQ OBJ  |
    #-----------------------------+
    my $inseq = Bio::SeqIO->new(-file   => "<$fasta_file",
				-format => 'fasta' ) 
	|| die "ERROR Can not open infile:\n $fasta_file\n";
 
    #-----------------------------+
    # OPEN GFF OUTPUT FILE        |
    #-----------------------------+
    open (GFFOUT, ">$gff_out") ||
	die "Can out open output file:\n$gff_out\n";

    #-----------------------------+
    # CHECK FOR GAPS IN THE OBJ   |
    #-----------------------------+
    while (my $seq = $inseq->next_seq()) {
	print "yup\n";
	#my $seq_string = $seq->seq();
	my $seq_len = $seq->length();
	#print "$seq_string\n";
	print "FULL SEQ LEN:\t$seq_len\n\n";

	# Increment across the seq string and see if this is 
	# the gap 
	for (my $i = 1; $i<$seq_len; $i++) {
	    my $seq_char = $seq->subseq($i,$i);
	    
	    #-----------------------------+
	    # DETERMINE IF THIS IS AS GAP |
	    # CHARACTER                   |
	    #-----------------------------+
	    if ( ($seq_char =~ $gap_char) || ($seq_char =~ "N") ) {
		#print "$i: $seq_char\n";
		$is_gap = 1;
	    }
	    else {
		$is_gap = 0;
	    }

	    #-----------------------------+
	    # DETERMINE START AND END OF  |
	    # GAPS                        |
	    #-----------------------------+
	    if ( ($is_gap) & (!$prev_gap) ) {
		# START OF A GAP
		$gap_start = $i;
	    } 
	    elsif ( (!$is_gap) & ($prev_gap) ) {
		# End of a gap
		$gap_end = $i;
		$gap_len = $gap_end - $gap_start;
		print "START:\t$gap_start\t";
		print "END:\t$gap_end\n";
		print "\tLEN: $gap_len\n";
		# If gap length is equal to or more then minimum 
		# write to gff file
		# This is labeled as gap
		if ($gap_len >= $min_gap_len) {
		    print "\tBig Enough\n";
		    
		    print GFFOUT 
			"$name_root\t".  # SeqName
			"gap\t".         # Source
			"gap\t".         # Feature (May need to make exon)
			"$gap_start\t".  # Start
			"$gap_end\t".    # End
			".\t".           # Score
			"+\t".           # Strand
			".\t".           # Frame
			"gap\n";         # Attribute

		}

	    } else {
		# Continuation of gap
	    }
	    

	    # Set the prev_gap for next round
	    if ($is_gap) {
		$prev_gap = 1;
	    }
	    else {
		$prev_gap = 0;
	    }
	    

	} # End of for $i
	
    } # End of while next_seq
    
    close GFFOUT;

    # TO DO: Convert the GFF output to Apollo game xml
    if (-e $gff_out) {
	apollo_convert($gff_out, "gff", $game_out, "game", 
		       $fasta_file, "NULL");
    }

} # End of for each file in the input folder

close LOG if $logfile;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    
    my $usage = "USAGE:\n".
	"  batch_findgaps.pl -i DirToProcess -o OutDir";

    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --len          # Min length to be considered gap\n".
	"  --gapchar      # Character indicating a gap\n".
	"  --logfile      # Path to file to use for logfile\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --test         # Run the program in test mode\n".
	"  --quiet        # Run program with minimal output\n";

    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
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
    $ApPath = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";

#    $ApPath = "/Applications/Apollo/bin/apollo";

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

=head1 NAME

batch_findgaps.pl - Annotate gaps in a fasta file

=head1 VERSION

This documentation refers to batch_findgaps version 1.0

=head1 SYNOPSIS

 Usage:
 batch_findgaps.pl -i DirToProcess -o OutDir

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

=item --m,mask

Single letter to mask with. Valid options are: [ N | n | X | x ]

=item --ext

The new outfile extension to use. Default value is .hard.fasta

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

The error messages that can be generated will be listed here.

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item *

RepeatMasker
(http://www.repeatmasker.org/)

=item *

Apollo (Genome Annotation Curation Tool)
http://www.fruitfly.org/annot/apollo/

=back

=head2 Required Perl Modules

=over

=item *

Getopt::Long

=back

=head1 BUGS AND LIMITATIONS

=head2 TO DO

=over 2


=item *

No current items on the to do list.

=back

=head2 Limitations

=over

=item *

No known majors limitations at this time.

=back

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/01/2007

UPDATED: 12/06/2007

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
