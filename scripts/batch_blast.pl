#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_blast.pl - Run blast on a set of fasta files.       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/23/2007                                       |
# UPDATED: 07/26/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of softmasked fasta files, this will   |
#  BLAST the files against a standard set of BLAST          | 
#  databases used for wheat genome annotation.              |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

batch_blast.pl - Run blast on a set of fasta files

=head1 VERSION

This documentation refers to batch_blast version 1.0

=head1 SYNOPSIS

 Usage:
 batch_blast.pl -i DirToProcess -o OutDir -d DbDir -c ConfigFile

=head1 DESCRIPTION

Given a directory of softmasked fasta files, this will
BLAST the files against a standard set of BLAST
databases used for wheat genome annotation.

All of the BLAST output files will be stored in a directory
name for the intput file in a subdirectory named blast.
(ie /home/myhome/infile/blast/infile_blastdb.blo).

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -d,--db-dir

Path of the directory containing the blast formatted databases.

=item -c, --config

Path to the batch_blast config file. This is a tab delimited text file
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

=item Col 2. Extension to use in blast output file. (ie. bln )

=item Col 3. Alignment output options (-m options from blast)

=item Col 4. Evalue threshold

=item Col 5. Database name

=item Col 6. Additional blast command switches

=back

An example config file would be:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=back

=head1 OPTIONS

=over 2

=item --blast-path

Full path to the NCBI blastall program. Default is blastall.

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

=item -verbose

Run the program with maximum output.

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

NCBI blastall
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST

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

Limited to NCBI BLAST

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

# PACKAGE LEVEL SCOPE

#Counters
my $i;                         # Number of blast files
my $count_db;                  # Number of databases to process
my $count_files;               # Number of files to process
my $count_proc;                # Number of processes
my $proc_num;                  # Number of the current process
my $file_num = 0;              # Number of the current file

# Files
my $file_config;               # Path to the Batch blast config file
my $file_to_blast;             # Path to the file to BLAST, the qry file
my $file_blast_out;            # Path to BLAST output file
my $logfile;                   # Path to a logfile to log error info
my $blast_path = "blastall";   # Path to the NCBI blastall binary
                               # by default will assume blastall works
                               # without providing full path

# Directories
my $dir_blast_out;             # Dir to hold the BLAST output for qry file 
my $dir_blast_db;
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $bac_out_dir;               # Dir for each sequnce being masked

my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $name_root;                 # Root name to be used for output etc

# ARRAYS
my @dbs;                       # Information for databases
                               # 2d Array filled from the config file

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required Arguments
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "d|db-dir=s" => \$dir_blast_db,
		    "c|config=s"   => \$file_config,
		    # Optional strings
		    "blast-path=s" => \$blast_path,
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

#//////////////////////
my $proc_num_max = 5;
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
if ( (!$indir) || (!$outdir) || (!$file_config) || (!$dir_blast_db) ) {
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
    print LOG "  batch_blast.pl\n";
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

unless ($dir_blast_db =~ /\/$/ ) {
    $dir_blast_db = $dir_blast_db."/";
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

$count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print "NUMBER OF FILES TO PROCESS: $count_files\n";


#-----------------------------+
# GET INFO FROM CONFIG FILE   |
#-----------------------------+
# Load to the 2d dbs array
open (CONFIG, $file_config) ||
    die "Can't open the config file:\n$file_config";
$i = 0;
my $line_num = 0;
while (<CONFIG>) {
    $line_num ++;

    unless (m/^\#/) {
	chomp;
	my @tmpary = split (/\t/);
	my $count_tmp = @tmpary;
	#print "COUNT: $count_tmp\n";
	#print "$_\n";
	#print "\t".$tmpary[0]."\n";
	if ($count_tmp == 6) {

	    #-----------------------------+
	    # TEST THE ENTRY              |
	    #-----------------------------+
	    #test_blast_db  ($prog, $db, $line, $dbdir)
	    test_blast_db( $tmpary[0], 
			   $tmpary[4], 
			   $line_num, 
			   $dir_blast_db );

	    $dbs[$i][0] = $tmpary[0];  # Blast Prog
	    $dbs[$i][1] = $tmpary[1];  # Outfile extension
	    $dbs[$i][2] = $tmpary[2];  # Alignment output
	    $dbs[$i][3] = $tmpary[3];  # E Value threshold
	    $dbs[$i][4] = $tmpary[4];  # DBNAME
	    $dbs[$i][5] = $tmpary[5];  # CMD SUFFIX
	    $i++;



	} # End of if count_tmp = 6
	else {
	    print "ERROR: Config file line number $lin_num\n";
	    print "       Only $lin_num variables were found\n" 
	}

    } # End of unless this is a comment line

} # End of while CONFIG


#-----------------------------+
# CHECK SANITY OF CONFIG FILE |
#-----------------------------+
# This will check for existence of blast database etc
print "Checking config file ...\n" unless $quiet;
for my $ind_db (@dbs) {
#    print "PROG:".$ind_db->[0]."\n";
#    print " EXT:".$ind_db->[1]."\n";
#    print "ALIG:".$ind_db->[2]."\n";
#    print "EVAL:".$ind_db->[3]."\n";
#    print "  DB:".$ind_db->[4]."\n";
#    print "SUFF:".$ind_db->[5]."\n";
#    print "\n\n";
}

$count_db = @dbs;
print "NUMBER OF DATABASES: $count_db\n";
$count_proc = $count_db * $count_files;
print "NUMBER OF PROCESSES: $count_proc\n";

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
# RUN BLAST FOR EACH FILE     |
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
    
    #-----------------------------+
    # Create parent dir if it     |
    # does not already exist      |
    #-----------------------------+
    my $dir_parent = $outdir.$name_root."/";
    unless (-e $dir_parent) {
	print "creating dir: $dir_parent\n";
	mkdir $dir_parent ||
	    die "Could not creat the output dir:\n$dir_parent\n";
    }


    #-----------------------------+
    # Create the dir to hold the  |
    # BLAST output                |
    #-----------------------------+
    $dir_blast_out = $outdir.$name_root."/blast/";
    unless (-e $dir_blast_out ) {
	print "Creating output dir\n: $dir_blast_out\n" if $verbose;
	mkdir $dir_blast_out ||
	    die "Could not create the output directory:\n$dir_blast_out";
    }
    
    $file_to_blast = $indir.$ind_file;


    #-----------------------------+
    # FOR EACH DB IN THE CONFIG   |
    # FILE                        |
    #-----------------------------+
    for my $ind_db (@dbs) {
	
	$proc_num++;

	# Temp exit for debug
	#if ($proc_num > $proc_num_max) { exit; }

	print "\n" unless $quiet;
	print "+-----------------------------------------------------------+\n"
	    if $verbose;
	print "| BLAST PROCESS: $proc_num of $count_proc \n" unless $quiet;
	print "+-----------------------------------------------------------+\n"
	    if $verbose;


	my $file_blast_out = $dir_blast_out.$name_root."_".$ind_db->[4]."."
	    .$ind_db->[1];

	my $blast_cmd = "$blast_path".
	    " -p ".$ind_db->[0].
	    " -i ".$file_to_blast.
	    " -d ".$dir_blast_db.$ind_db->[4].
	    " -e ".$ind_db->[3].
	    " -o ".$file_blast_out.
	    " -m ".$ind_db->[2].
	    " ".$ind_db->[5];
	
	# PRINT VERBOSE OUTPUT
	if ($verbose) {
	    print "\tFILE: $name_root\n";
	    print "\t  DB: ".$ind_db->[4]."\n";
	    print "\t CMD: $blast_cmd\n";
	}

	unless ($test) {
	    system ($blast_cmd);
	}



    } # End of for each ind_db



} # End of for each file in the input folder

close LOG if $logfile;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub test_blast_db {
# Check for the existence of a BLAST database

    # prog ----- blast program to use
    # $db ------ blast database name
    # $line ---- line number from the config file
    # $dbdir --- path to the blast database directory
    # $db_path - path to the expected blast database

    my ($prog, $db, $line, $dbdir) = @_;
    my $db_path;
   
    #-----------------------------+
    # Nucleotide BLAST            |
    #-----------------------------+
    if ( ($prog =~ "blastn") || ($prog =~ "megblast") ||
	 ($prog =~ "tblastx") || ($prog =~ "tblastn")  ) {
	$db_path = $dbdir.$db.".nin";
    }
    #-----------------------------+
    # Protein BLAST               |
    #-----------------------------+
    elsif ( ($prog =~ "blastx") || ($prog =~ "blastp") ||
	    ($prog =~ "psi-blast") || ($prog =~ "phi-blast") ) {
	# Protein BLAST
	$db_path = $dbdir.$db.".pin";
	
    }
    #-----------------------------+
    # Unrecognized BLAST          |
    #-----------------------------+
    else {
	print "\a";
	print "ERROR: Config file line $line\n";
	print "The blast program $prog is not recognized by batch_blast.pl\n";
	exit;
    }

    #-----------------------------+
    # CHECK FOR DB FILE EXISTNECE |
    #-----------------------------+
    unless (-e $db_path) {
	print "\a";
	print "ERROR: Config file line $line\n";
	print "The database file for $db could not be found at:\n".
	    "$db_path\n";
    } else {
	print "The db file is okay at: $db_path\n" if $verbose;
    }

}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    
    my $usage = "USAGE:\n".
	"  batch_blast.pl -i DirToProcess -o OutDir";

    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"  --config       # Path to the batch_blast config file\n".
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

STARTED: 07/23/2007

UPDATED: 07/24/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/23/2007
# - Program started
# - Imported functions from batch_mask.pl as well as
#   Blast_PAWS.pl
# - Added use of a config file to define attributes of
#   BLAST jobs for each database to BLAST against
#   where lines starting with # are ignored
#
# 07/24/2007
# - Adding test_blast_db to test the existence of blast
#   databases
