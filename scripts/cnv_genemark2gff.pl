#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_repmask2gff.pl - Convert repeatmasker output to gff   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_sourceforge.net                   |
# STARTED: 10/30/2007                                       |
# UPDATED: 10/30/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Runs the RepeatMasker program for a set of input         |
#  FASTA files against a set of repeat library files &      |
#  then converts the repeat masker *.out file into the      |
#  GFF format and then to the game XML format for           |
#  visualization by the Apollo genome anotation program.    |
#                                                           |
# USAGE:                                                    |
#  cnv_repmaske2gff.pl -i infile.out -o outfile.gff         |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;
print "\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
use Getopt::Long;
use Bio::Tools::Genemark;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 259 $ =~ /(\d+)/;
#//////////////////////
my $file_num_max = 1;
#\\\\\\\\\\\\\\\\\\\\\\

#-----------------------------------------------------------+
# VARIABLE SCOPE                                            |
#-----------------------------------------------------------+

# VARS WITH DEFAULT VALUES
my $engine = "crossmatch";

# GENERAL PROGRAM VARS
my @inline;                    # Parse of an input line
my $msg;                       # Message printed to the log file

# DIR PATHS
my $infile;                    # Input path of *.out file to convert
my $outfile;                   # Output gff file path
my $prefix;                    # Prefix name for gff output file

my $bac_out_dir;               # Dir for each sequence being masked
my $bac_rep_out_dir;           # Dir to hold BAC repeat masker output
                               # $bac_out_dir/rm
my $gff_out_dir;               # Dir for the gff output
                               # $bac_out_dir/gff

# FILE PATHS
my $logfile;                   # Path to a logfile to log error info
my $rm_path;                   # Full path to the repeatmasker binary
my $ap_path;                   # Full path to the apollo program
my $rep_db_path;               # Path to an indivual repeat database
my $file_to_mask;              # The fasta file to be masked
my $repmask_outfile;           # Repeat masked outfile
my $repmask_catfile;           # Concatenated repeat masked outfile
my $config_file;               # Full path to the configuration file
                               # This includes the db names and paths
                               # of the fasta files to use for masking
my $repmask_cat_cp;            # Path to the copy of the RepMask 



my $src_prog = "GeneMarkHMM";        # Source program and matrix
my $src_seq = "unknown";

my $name_root;                 # Root name to be used for output etc
my $repmask_log;
my $repmask_log_cp;

my $repmask_tbl_file;
my $repmask_masked_file;


my $gff_alldb_out;
my $gff_el_out;
my $xml_alldb_out;
my $xml_el_out;

# FINAL FILE LOCATIONS
my $repmask_masked_cp;        # Masked fasta file
my $repmask_local_cp;         # Copy of masked in fasta file
my $repmask_tbl_cp;           # Repmask Table
my $repmask_el_cp;            # Individual 
my $repmask_xml_el_cp;
my $repmask_out_cp;           # Repmask out file copy

# REPEAT DB VARS
my $rep_db_name;               # Name of the repeat database
my $ind_lib;                   # Vars for an individual library

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;
my $debug = 0;                 # Run the program in debug mode 

# PROGRAM COMMAND STRINGS
my $cmd_repmask;               # Command to run RepeatMasker
my $cmd_make_gff_db;           # Make the gff file for an individual database
my $cmd_make_gff_el;           # Appears to not be used

# COUNTERS AND INDEX VARS
my $num_proc = 1;              # Number of processors to use
my $i;                         # Used to index the config file
my $proc_num = 0;              # Counter for processes
my $file_num = 0;              # Counter for fasta files being processed

# ARRAYS
my @mask_libs = ();

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # REQUIRED ARGUMENTS
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    "prefix"       => \$prefix,
		    # ADDITIONAL OPTIONS
		    "src-prog=s"   => \$src_prog, 
		    "src-seq=s"    => \$src_seq,
		    "rm-path=s"    => \$rm_path,
		    "ap-path=s",   => \$ap_path,
		    "logfile=s"    => \$logfile,
		    "p|num-proc=s" => \$num_proc,
		    "engine=s"     => \$engine,
		    # BOOLEANS
		    "apollo"       => \$apollo,
		    "verbose"      => \$verbose,
		    "debug"        => \$debug,
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
	"Version: $VERSION\n\n";
    exit;
}

# Show full help when required options
# are not present
if ( (!$infile) || (!$outfile) ) {
    print "\a";
    print "ERROR: An input file was not specified with the -i flag\n"
	if !$infile;
    print "ERROR: An output file was not specified with the -o flag\n"
	if !$outfile;
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
    print LOG "  batch_mask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
}

genemark_to_gff ($infile, $outfile, $src_seq, $src_prog );

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
	"  batch_mask.pl -i DirToProcess -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"  --config       # Path to database list config file\n".
	"\n".
	"ADDITIONAL OPTIONS:\n".
	"  --rm-path      # Full path to repeatmasker binary\n".
	"  --engine       # The repeatmasker engine to use:\n".
	"                 # [crossmatch|wublast|decypher]\n".
	"                 # default is to use crossmatch\n".
	"  --num-proc     # Number of processors to use for RepeatMasker\n".
	"                 # default is one.\n".
	"  --apollo       # Convert output to game.xml using apollo\n".
	"                 # default is not to use apollo\n".
	"  --quiet        # Run program with minimal output\n".
	"  --test         # Run the program in test mode\n".
	"  --logfile      # Path to file to use for logfile\n".
	"\n".
	"ADDITIONAL INFORMATION\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n";

	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}

sub genemark_to_gff {
    
    my ($gm_in_path, $gff_out_path, $gm_src_seq, $gm_src_prog) = @_;

    # OPEN THE GENEMARK INFILE
    my $gm_obj = Bio::Tools::Genemark->new(-file => $gm_in_path);

    # OPEN THE GFF OUTFILE
    my $rna_count = 0;
    while(my $gene = $gm_obj->next_prediction()) {
       
	$rna_count++;
	#$result = sprintf("%08d", $number);
	my $rna_id = sprintf("%04d", $rna_count);

	my @exon_ary = $gene->exons();
	my $num_exon = @exon_ary;

	#print "START\tEND\tORIENT";
	for my $ind_gene (@exon_ary) {
	    my $start = $ind_gene->start;
	    my $end = $ind_gene->end;
	    my $strand = $ind_gene->strand;
	    if ($strand =~ '1') {
		$strand = "+"; 
	    }
	    elsif ($strand =~ '-1') {
		$strand = "-";
	    }
	    else {
		$strand = ".";
	    }
	    
	    # Pseudo gff
	    print $gm_src_seq."\t".   # seq name
		$gm_src_prog."\t".    # source
		"exon\t".             # feature
		$start."\t".          # start
		$end."\t".            # end
		".\t".                # score
		$strand."\t".         # strand
		".\t".                # frame
		"RNA$rna_id\n";                 # attribute
	}

    }

    $gm_obj->close();

}

=head1 NAME

batch_mask.pl - Run RepeatMasker and parse results to a gff format file. 

=head1 VERSION

This documentation refers to program version $Rev: 259 $

=head1 SYNOPSIS

 USAGE:
    batch_mask.pl -i DirToProcess -o OutDir -c ConfigFile

=head1 DESCRIPTION

Runs the RepeatMasker program for a set of input
FASTA files against a set of repeat library files &
then converts the repeat masker *.out file into the
GFF format and then to the game XML format for
visualization by the Apollo genome anotation program.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Configuration file that lists the database names and paths of the
fasta files to use as masking databases.

=back

=head2 Additional Options

=over 2

=item -p,--num-proc

The number of processors to use for RepeatMasker. Default is one.

=item --engine

The repeatmasker engine to use: [crossmatch|wublast|decypher].
The default is to use the crossmatch engine.

=item --apollo

Use the apollo program to convert the file from gff to game xml.
The default is not to use apollo.

=item --rm-path

The full path to the RepeatMasker binary.

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands.

=back

=head2 Additional Program Information

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

=back

=head1 DIAGNOSTICS

The error messages that can be generated will be listed here.

=head1 CONFIGURATION AND ENVIRONMENT

The major configuration file for this program is the list of
datbases indicated by the -c flag.

=head2 Databases Config File

This file is a tab delimited text file. Line beginning with # are ignored.

B<EXAMPLE>

  #-------------------------------------------------------------
  # DBNAME       DB_PATH
  #-------------------------------------------------------------
  TREP_9         /db/repeats/Trep9.nr.fasta
  TIGR_Trit      /db/repeats/TIGR_Triticum_GSS_Repeats.v2.fasta
  # END

The columns above represent the following 

=over 2

=item Col. 1

The name of the repeat library
This will be used to name the output files from the analysis
and to name the data tracks that will be used by Apollo.

=item Col. 2

The path to the fasta format file containing the repeats.

=back

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

=item *

Make it a variable to possible to put the gff output (1) all in positive 
strand, (2) all in negative strand, (3) alignment to positive or
negative strand, (4) cumulative in both positive and negative strand.
Current behavior is to do number 1 above.

=back

=head2 Limitations

=over

=item *

This program has been tested with RepeatMasker v  3.1.6

=back

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 04/10/2006

UPDATED: 09/11/2007

VERSION: $Rev: 259 $

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
# 07/13/2007
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
# 09/07/2007
# - Moving POD documentation to the end of the program
# - Changing to use a config file instead of internal 2-d array
# - Modified all variable names to lowercase
# - Added use strict
#
# 09/11/2007
# - Getting rid of LOG, all printing to STDERR
# - Dropped attempts to move *.tbl file.
