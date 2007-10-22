#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_dirsplit.pl - split dir of fasta files to n dirs    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/24/2007                                       |
# UPDATED: 07/24/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a fasta files, moves or copies the files into n    |
#  separate dirs where n is a variable passed at the        |
#  command line.                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

fasta_dirsplit.pl - Split dir of fasta files into n dirs

=head1 VERSION

This documentation refers to fasta_dirsplit version 1.0

=head1 SYNOPSIS

 Usage:
 fasta_dirsplit.pl -i DirToProcess -o OutDir -b BaseName -n NumDirs

=head1 DESCRIPTION

Spilts a directory containing multiple fasta files into n directories.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -n,--num-dir

Number of dirs to split the parent dir into.

=item -b, --base-name

Name to use a base name for creating the subdirectory

=item 

=back

=head1 OPTIONS

=over 2

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
#my $out_ext = ".hard.fasta";  # Outfile extension
#my $mask_char = "N";          # Character to mask with

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode

# PACKAGE LEVEL SCOPE
my $num_dir;                   # Number of dirs to split into
my $file_to_move;              # Path to the file to mask
my $file_new_loc;              # New location of the file to moave
#my $hard_mask_out;             # Path to the hardmasked file output
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $base_name;                 # Base name to be used for output etc

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "n|num-dir=s"   => \$num_dir,
		    "b|base-name=s" => \$base_name,
		    # Optional strings
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);



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

unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# CREATE OUT DIR SUBDIR       |
# IF IT DOES NOT EXIST        |
#-----------------------------+
my $max_num = $num_dir + 1;
for (my $i=1; $i<$max_num; $i++) {
    my $newsubdir = $outdir.$base_name.$i;
    
    unless (-e $newsubdir) {
	print "Making dir $newsubdir\n" if $verbose;
	mkdir $newsubdir ||
	    die "Could not create the output directory:\n$newsubdir\n"
    }
    else {
	print "WARNING: The directory already exists:\n\t$newsubdir\n";
	print "\tOutput will be appended to this direcotry\n";
    }
    
}


# Temp exit
#exit;

#-----------------------------+
# RUN REPEAT MAKSER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# RepLibs ARRAY               |
#-----------------------------+

$i = 1;
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

    $file_to_move = $indir.$ind_file;
    $file_new_loc = $outdir.$base_name.$i."/".$ind_file;
    
    print "\nCopying:\n".
	"\t$file_to_move TO\n".
	"\t$file_new_loc\n" if $verbose;
    my $errmsg = "Can not move file:\n\t".
	"FROM:$file_to_move\n\t".
	"  TO:$file_new_loc\n";
    copy ($file_to_move, $file_new_loc) ||
	die "$errmsg\n";
    
    if ($i == $num_dir) {
	$i = 1;
    }
    else {
	$i++;
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
	" fasta_dirsplit.pl -i DirToProcess -o OutDir -b BaseName".
	" -n NumDirs\n";

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

STARTED: 07/24/2007

UPDATED: 07/25/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/24/2007
# - Program started
# - Basic outline with help, usage and man working
# - Reads dir of fasta files and moves names to array
