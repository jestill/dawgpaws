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

=head1 HISTORY

STARTED: 07/19/2007

UPDATED: 07/19/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/19/2007
# - Program started
# - Imported functions from batch_mask.pl as well as
#   HardMask.pl
