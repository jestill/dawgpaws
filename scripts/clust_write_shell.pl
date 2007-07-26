#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# clust_write_shell.pl - Write shell scripts for r cluster  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/26/2007                                       |
# UPDATED: 07/26/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given information from the command line, write the shell |
#  scripts required to run jobs on the rcluster.            |
#                                                           |
# LICENSE                                                   |
#  GNU LESSER GENERAL PUBLIC LICENSE                        |
#  http://www.gnu.org/licenses/lgpl.html                    |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

clust_write_shell.pl - Write shell scripts for the r cluster

=head1 VERSION

This documentation refers to fasta_dirsplit version 1.0

=head1 SYNOPSIS

 Usage:
 clust_write_shell.pl -p program -o outdir -b job -n 16
  -p   # program to write the shell script for
  -o   # output dirctory to 
  -b   # base name for the shell scripts to write
  -n   # number of shell scripts to write
       # this will also be used to reference dirs

=head1 DESCRIPTION

Spilts a directory containing multiple fasta files into n directories.

=head1 REQUIRED ARGUMENTS

=over 2

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

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $submit_job = 0;            # Submit the shell script that is created
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
my $prog_name;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "p|program=s"   => \$prog_name,
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
if ( (!$outdir) || (!$prog_name) || (!$num_dir) || (!$base_name) ) {
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
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
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
# WRITE SHELL SCRIPTS         |
#-----------------------------+

my $max_num = $num_dir + 1;
for (my $i=1; $i<$max_num; $i++) {

    # open the shell script
    my $shell_path = $outdir.$base_name."_shell".$i.".sh";
    
    print "Writing to $shell_path\n";
    #open ()
    if ($prog_name =~ "batch_blast") {
	open (SHOUT, ">".$shell_path);
	
	print SHOUT "/home/jlblab/jestill/dawg-paws/scripts/batch_blast.pl".
	    " -i /scratch/jestill/wheat_in/$base_name$i/".
	    " -o /scratch/jestill/wheat_out/".
	    " -c /home/jlblab/jestill/scripts/dawg-paws/batch_blast_full.jcfg".
	    " -d /db/jlblab/".
	    " --logfile /home/jlblab/jestill/$base_name$i.log";
    } 
    else {
	print "The program name is not recognized: $prog_name\n";
	exit;
    }
    
}



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

STARTED: 07/26/2007

UPDATED: 07/26/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/24/2007
# - Program started
# - Basic outline with help, usage and man working
# - Reads dir of fasta files and moves names to array
