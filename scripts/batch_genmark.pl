#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_gemark.pl - Run Genmark gene prediction program     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com                         |
# STARTED: 11/09/2007                                       |
# UPDATED: 11/09/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the genmark gene prediction program in batch mode.   |
#  Runs genmark as well as converts output to gff format.   |
#  Requires a config file to specify libraries to use.      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use File::Copy;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $file_config;               # Path to the configuration file
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file

my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc

# GENMARK PATHS
# genmakr_dir is the dir that contains the genmark binaries
my $genmark_dir = $ENV{GM_BIN_DIR};
# lib_dir is the dir that contains the genmark matrix librarires
my $lib_dir = $ENV{GM_LIB_DIR};

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;

# COUNTERS
my $num_proc = 1;              # Number of processes

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "c|config=s"    => \$file_config,
		    # Optional strings
		    "genmark-dir=s" => \$genmark_dir,
		    "lib-dir=s"     => \$lib_dir,
		    "logfile=s"     => \$logfile,
		    # Booleans
		    "apollo"        => \$apollo,
		    "verbose"       => \$verbose,
		    "test"          => \$test,
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,
		    "q|quiet"       => \$quiet,);

my $proc_num = 0;

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
    print "\batch_genmark.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

# Show full help when required options
# are not present
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print_help("full");
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


for my $ind_file (@fasta_files)
{

    $proc_num++;
    $file_num++;
    #if ($file_num == $file_num_max){exit;}

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	# This is hard masked fasta files
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE GENMARK OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $genmark_out_dir = $outdir.$name_root."/genemark/";
    unless (-e $genmark_out_dir) {
	mkdir $genmark_out_dir ||
	    die "Could not create genscan out dir:\n$genmark_dir\n";
    }
    
    my $infile_path = $indir.$ind_file;
    my $out_path = $genmark_out_dir.$name_root.".genmark.out";
    my $gff_path = $genmark_out_dir.$name_root.".genmark.gff";

    # To be more generalizable the following vars for each matrix
    # should be read in from a config file, but I am cutting 
    # corners here and hard coding this
    my $gm_os_cmd = $genmark_dir."gmhmme3 -m ".$lib_dir."o_sativa.mod".
	" -o $out_path $infile_path";
    my $gm_zm_cmd = $genmark_dir."gmhmme2 -m ".$lib_dir."corn.mtx".
	" -o $out_path $infile_path";
    my $gm_ta_cmd = $genmark_dir."gmhmme2 -m ".$lib_dir."wheat.mtx".
	" -o $out_path $infile_path";
    my $gm_hv_cmd = $genmark_dir."gmhmme2 -m ".$lib_dir."barley.mtx".
	" -o $out_path $infile_path";

    print "=======================================\n" if $verbose;
    print "Running Genmark for $name_root\n" if $verbose;
    print " File $file_num of $num_files\n" if $verbose;
    print "=======================================\n" if $verbose;

    print "\n$gm_os_cmd\n" if $verbose;
    system($gm_os_cmd) unless $test;

    print "\n$gm_zm_cmd\n" if $verbose;
    system($gm_zm_cmd) unless $test;

    print "\n$gm_ta_cmd\n" if $verbose;
    system($gm_ta_cmd) unless $test;

    print "\n$gm_hv_cmd\n" if $verbose;
    system($gm_hv_cmd) unless $test;


#    print "Converting Genscan output\n";
#    if (-e $out_path) {
#	genscan_2_gff($out_path, $gff_path, $name_root) unless $test;
#    }
#    else {
#	print "ERROR: Could not find genscan output at:\n$out_path\n"
#    }

} # End of for each file in the input folder

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub genmark_2_gff 
{

} #End of genscan_2_gff subfunction


sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;
    
    my $usage = "USAGE:\n".
	"  batch_genscan.pl -i DirToProcess -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --gencan-path  # Full path to the genscan binary\n".
	"  --lib-path     # Full path to the prediction library:\n".
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


=head1 NAME

batch_genmark.pl - Run GenMark.hmm and parse results to a gff format file. 

=head1 VERSION

This documentation refers to batch_mask version 1.0

=head1 SYNOPSIS

 Usage:
 batch_genscan.pl -i DirToProcess -o OutDir

=head1 DESCRIPTION

Run the genscan gene prediction program in batch mode.
This will run genscan as well as convert the output to gff format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --genscan-path

The full path to the genscan binary.

=item --lib-path

The full path to the library file.

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

Genscan

=back

=head2 Required Perl Modules

=over

=item *

File::Copy

=item *

Getopt::Long

=back

=head1 BUGS AND LIMITATIONS

=head2 TO DO

=over 2

=item *

Make the results compatable for an upload to a chado
database.

=back

=head2 Limitations

=over

=item *

Currently no known limitations.

=back

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 11/09/2007

UPDATED: 11/09/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 11/09/2007
# - Program started to run GenMark.hmm.euk in batch mode
#   on the local machine.
# 
