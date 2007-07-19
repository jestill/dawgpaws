#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_hardmask.pl - Hardmask a batch of softmasked files  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/19/2007                                       |
# UPDATED: 07/19/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of softmasked fasta files, this will   |
#  hardmask the files by replacing the lowecase letters     |
#  with an uppercase letter set by the user.                |
#                                                           |
# USAGE:                                                    |
#  batch_hardmask.pl -i InDir -o OutDir -m X                |
#                                                           |
# LICENSE                                                   |
#  GNU LESSER GENERAL PUBLIC LICENSE                        |
#  http://www.gnu.org/licenses/lgpl.html                    |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

batch_hardmask.pl - Hardmask a batch of softmasked fasta files. 

=head1 VERSION

This documentation refers to batch_hardmask version 1.0

=head1 SYNOPSIS

 Usage:
 batch_hardmask.pl -i DirToProcess -o OutDir

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

=cut

print "\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use File::Copy;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $ver = "1.0";

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


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "m|mask"       => \$mask_char,
		    "ext"          => \$out_ext,
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
my $file_num_max = 5;
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
    print LOG "  batch_hardmask.pl\n";
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

print "Creating output dir ...\n" unless $quiet;
unless (-e $outdir) {
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# RUN REPEAT MAKSER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# RepLibs ARRAY               |
#-----------------------------+

for my $ind_file (@fasta_files)
{
    
    $file_num++;


    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
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

    $file_to_mask = $indir.$ind_file;
    $hard_mask_out = $outdir.$name_root.$out_ext;
    
    print "Converting:\n".
	"\t$file_to_mask TO\n".
	"\t$hard_mask_out\n" unless $quiet;


    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (IN, "<".$file_to_mask) ||
	die "Can not open input file:\n$file_to_mask\n";
    
    open (OUT, ">".$hard_mask_out) ||
	die "Can not open output file:\n$hard_mask_out\n";
    
    #-----------------------------+
    # HARD MASK FILE              |
    #-----------------------------+
    # The tr regexp does not appear to accept variables
    # therefore I have to write this a bit convoluted with
    # if then statements for acceptable MaskCharacters
    while (<IN>)
    {
	unless (m/^>/)        # Do not mask header lines 
	{
	    # Mask with the selected character
	    if ($mask_char =~ "N"){
		tr/[a-z]/N/;
	    } elsif ($mask_char =~ "X"){
		tr/[a-z]/X/;
	    } elsif ($mask_char =~ "x"){
		tr/[a-z]/x/;
	    } elsif ($mask_char =~ "n"){
		tr/[a-z]/n/;
	    } else {
		$msg = "\aERROR: A valid mask character was not selected\n";
	    }# End of select mask character
		
		# Print masked string to the outfile
		print OUT $_;
	    
	} else {
	    print OUT $_;
	    $NumRecs++;       # For headers increment NumRecs
	}
    } # End of while IN
    
    #-----------------------------+
    # CLOSE FILES                 |
    #-----------------------------+
    close IN;
    close OUT;
    
#//////////////////////////////////
# MAY WANT TO LEAVE THE FOLLOWING
# TO MAKE A COPY IN THE GENERAL
# BAC DIRECTORY
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#
##    #-----------------------------+
#    # MAKE OUTPUT DIR             |
##    #-----------------------------+
#    # the $bac_out_dir is the summary directory that
#    # contains all of the data related to a seq
#    $bac_out_dir = $outdir.$name_root."/";
#    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 


#    # TEMP EXIT FOR DEBUG, WIll JUST RUN FIRST FILE TO BE MASKED
#    if ($file_num > $file_num_max ) {
#	print "\nDebug run finished\n\n";
#	exit;
#    }


} # End of for each file in the input folder

close LOG if $logfile;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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

    print "Converting output to Apollo game.xml format\n"
	unless $quiet;

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

    #$ApPath = "/home/jestill/Apps/Apollo/";
    $ApPath = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
#    $ApPath = "/home/jestill/Apps/Apollo/";
#    $ApCmd = $ApPath."Apollo -i ".$InForm." -f ".$InFile.
#	" -o ".$OutForm." -w ".$OutFile;

    $ApCmd = $ApPath." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" ) {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb") {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }
    
    # Do the apollo command
    system ( $ApCmd );

}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    
    my $usage = "USAGE:\n".
	"  batch_hardmask.pl -i DirToProcess -o OutDir";

    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --mask         # Character to mask with [N|n|X|x]\n".
	"                 # default is N\n".
	"  --ext          # Extension to add to new files.\n".
	"                 # default is .hard.fasta.\n".
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

sub rmout_to_gff {

# Subfunction to convert repeatmasker out file
# to gff format for apollo
    
    # $rm_file = path to the repeat masker file
    # $gff_file = path to the gff output file
    # $pre is the way that the output file will
    # be made, it should be either > or >>
    # > for overwrite
    # >> for concatenate
    # the default will be to concatenate when the prefix
    # variable can not be undertood
    my ( $rm_file, $gff_file, $pre ) = @_;
    my $IN;
    my $OUT;
    my $strand;

    #-----------------------------------------------------------+
    # REPEATMASKER OUT FILE CONTAINS
    #-----------------------------------------------------------+
    # 0 Smith-Waterman score of the match
    # 1 % substitutions in matching region compared to the consensus
    # 2 % of bases opposite a gap in the query sequence (deleted bp)
    # 3 % of bases opposite a gap in the repeat consensus (inserted bp)
    # 4 name of query sequence
    # 5 starting position of match in query sequence
    # 6 ending position of match in query sequence
    # 7 no. of bases in query sequence past the ending position of match
    #   in parenthesis
    # 8 match is with the Complement of the consensus sequence in the database
    #   C
    # 9 name of the matching interspersed repeat
    #10 class of the repeat, in this case a DNA transposon 
    #11 no. of bases in (complement of) the repeat consensus sequence 
    #   in parenthesis
    #12 starting position of match in database sequence 
    #   (using top-strand numbering)
    #13 ending position of match in database sequence

    open ( RM_IN, $rm_file ) ||
	die "Could not open the RM infile\n";

#    open ( RM_OUT, ">".$gff_file) ||
    open ( RM_OUT, $pre.$gff_file) ||
	die "Could not open the GFF outfile\n";
    
    while (<RM_IN>) {
	chomp;
	my @rmout = split;

	my $line_len = @rmout;

	my $cur_strand = $rmout[8];

	if ($cur_strand =~ "C" ) {
	    $strand = "-";
	}
	else {
	    $strand = "+";
	}

	#-----------------------------+
	# The following uses the TREP |   
	# hits as names this is done  |
	# by calling the attributes   |
	# exons..a strange kluge      |
	#-----------------------------+
	# This places the RepMask output in the same frame as 
	# the hit was found in the database. 
	#print RM_OUT "$rmout[4]\t";      # qry sequence name
	#print RM_OUT "repeatmasker:trep9\t";   # software used
	#print RM_OUT "exon\t";  # attribute name
	#print RM_OUT "$rmout[5]\t";      # start
	#print RM_OUT "$rmout[6]\t";      # stop
	#print RM_OUT "$rmout[0]\t";      # smith waterman score"
	#print RM_OUT "$strand\t";              # strand
	#print RM_OUT ".\t";              # frame
	#print RM_OUT "$rmout[9]";        # attribute
	#print RM_OUT "\n";

	#-----------------------------+
	# THE FOLLOWING MAKES THE HITS|
	# CUMULATIVE IN BOTH STRANDS  |
	# OR JUST THE POSITVE STRAND  |
	#-----------------------------+
	print RM_OUT "$rmout[4]\t";      # qry sequence name
	print RM_OUT "repeatmasker:trep9\t";   # software used
	print RM_OUT "exon\t";  # attribute name
	print RM_OUT "$rmout[5]\t";      # start
	print RM_OUT "$rmout[6]\t";      # stop
	print RM_OUT "$rmout[0]\t";      # smith waterman score"
	print RM_OUT "+\t";              # Postive strand
	print RM_OUT ".\t";              # frame
	print RM_OUT "$rmout[9]";        # attribute
	print RM_OUT "\n";

	#print RM_OUT "$rmout[4]\t";      # qry sequence name
	#print RM_OUT "repeatmasker:trep9\t";   # software used
	#print RM_OUT "exon\t";  # attribute name
	#print RM_OUT "$rmout[5]\t";      # start
	#print RM_OUT "$rmout[6]\t";      # stop
	#print RM_OUT "$rmout[0]\t";      # smith waterman score"
	#print RM_OUT "-\t";                # Negative strand
	#print RM_OUT ".\t";              # frame
	#print RM_OUT "$rmout[9]";        # attribute
	#print RM_OUT "\n";

    }

    return;
    
}


=head1 HISTORY

STARTED: 04/10/2006

UPDATED: 07/18/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 4/10/2006
# - Program started
# - Parsing of RepeatMasker *.out file to an Apollo 
#   compatible *.GFF file done. Working copy for both
#   a single track for each repeat class, and a data
#   track with all of the results from a Repeat Library
#   grouped toether.
# - Two dimensional array containing the repeat library
#   database informaiton.
#
# 4/11/2006
# - Additional code comments and reformat
# - Code to cycle through a set of FASTA files and store
#   the output data in a separate named folder for each
#   of the FASTA flies.
# - Added the -fixed to the command to run RepeatMasker
#   to make sure that the text format is the same for all
# - Since Apollo does not allow additional GFF files to be
#   added later, I added code to create a GFF file that has
#   the outcome for all of the RepeatMask runs in one 
#   database.
#
# 4/12/2006
# - Added the AllRepeats to the dataset to make sure that 
#   that the repeat characterization is cumulative
# - ConvertGFF2Chado subfunction
# - Adding Smith Waterman score from RepeatMasker to the
#   GFF output file. This should be the first column in the 
#   *.out file and the 6th column in the GFF file.
#
# 4/16/2005
# - Added array for list of fasta files to process and testing
#   with small set of sequences.
# - The fasta file was manually edited to use just the GB
#   ID in the FASTA header. This could be done using the 
#   the seq object PERL module from bioperl.
# - The fasta files should be read from an input directory
#   and then the output should be copied to an output dir.
#
# 4/17/2006
# - Working out variables names
#
# 4/25/2006
# - Adding ability to get the input set of FASTA files
#   from a directory 
#
# 4/26/2006
# - Working out the copy of the relevant output files to
#   a central directory for the FASTA file (ie. all 
#   programatic output from the FASTA file goes to 
#   a dir named for the file.) 
# - Added use of File::Copy module for cp commands
#
# 5/24/2006
# - Added process log to the program. This will allow
#   me to launch the program at the lab and then
#   monitor the process at home.
#
# 5/31/2006
# - Made local copy dir a variable that can be set at the top
#
# 09/26/2006
# - A few changes made to the base code to make this
#   work on the altix. Using databases at
#   /scratch/jestill/repmask/
#
# 09/28/2006
# - Changes make to get this to work on the altix
# - Changed the format of the repeat databases list
#   This is currently a two-d array
#
#-----------------------------+
# 07/13/2007 - VERSION 1.0    |
#-----------------------------+
# - Added POD documentation.
# - Renamed to automask.pl
# - Made this the official 1.0 release
# - Adding command line variables
# - Added print_help subfunction
# - Modified ApolloConvert subfunction to 
#   apollo_convert
# - Getting rid of the die commands and replacing
#   them with the write to logfile commands. This
#   should prvent the program from dying when
#   working on really large batch files.
# - Added cmd line options for:
#   - Number of processors
#   - Engine to use in repeatmasker
#   - apollo variable at cmd line
#     initially will be used as boolean but can be
#     used to pass apollo version, path and
#     desired output (ie. chado, game.xml etc)
#
# 07/15/2007
# - Added the ability to show help when the
#   required options $indir and $outdir are not
#   present.
#
# 07/16/2007
# - Added error message when no fasta files were
#   found in the input directory.
# - Added the ability to append a slash to the end
#   of the $indir and $outdir variables if they
#   are not already present
#
# 07/17/2007
# - Added the ability to get the search_name directly
#   from the fasta file header
# - Getting root name for the output files and folders
#   by stripping off the fasta extension with regexp
# - Changed output dir to just the root name
# - Adding rmout_to_gff subfunction to convert the repeat
#   masker output to gff format. This will get rid
#   of any dependence on grep and awk and will provide
#   a more stable program
# 
# 07/18/2007
# - rmout_to_gff working now ... strangely had to make
#   the type exon to draw properly in apollo
# - Got rid of previous code that used grep and awk
#   to do this  conversion
# - Added command line variable to specify the full
#   path of the RepeatMasker program. This should
#   help this to run machines without having to 
#   make any changes to the user's environment.
#   (ie. don't have to add RepeatMaker to user's path)
# - Renamed program again to batch_mask.pl
# 
