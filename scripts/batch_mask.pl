#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_mask.pl - Run RepeatMasker in batch mode            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_sourceforge.net                   |
# STARTED: 04/10/2006                                       |
# UPDATED: 09/11/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Runs the RepeatMasker program for a set of input         |
#  FASTA files against a set of repeat library files &      |
#  then converts the repeat masker *.out file into the      |
#  GFF format and then to the game XML format for           |
#  visualization by the Apollo genome anotation program.    |
#                                                           |
# USAGE:                                                    |
#  batch_mask.pl -i InDir -o OutDir -c ConfigFile.txt       |
#                                                           |
#  batch_mask.pl --man                                      |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;
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
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
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


my $search_name;               # Name searched for in grep command
my $name_root;                 # Root name to be used for output etc
my $repmask_log;
my $repmask_log_cp;

my $repmask_tbl_file;
my $repmask_masked_file;


my $gff_alldb_out;
my $gff_el_out;
my $xml_alldb_out;
my $xml_el_out;
#my $repmask_cat_cp;

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
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "c|config=s",  => \$config_file,
		    # ADDITIONAL OPTIONS
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

my $bac_parent_dir = $outdir;  

my ( @rep_libs );


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
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print "ERROR: An input directory was not specified with the -i flag\n"
	if !$indir;
    print "ERROR: An output directory was not specified with the -o flag\n"
	if !$outdir;
    print "ERROR: A config file was not specified with the -c flag\n"
	if !$config_file;
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

#-----------------------------+
# REPEAT LIBRARY INFORMATION  |
#-----------------------------+
open (CONFIGFILE, "<$config_file") 
    || die "Could not open the config file:\n$config_file\n";
    

$i = 0;

my $config_line_num = 0;
print "Parsing config file...\n" if $verbose;
while (<CONFIGFILE>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
	# Split input by tab 
	my @in_line = split(/\t/, $_);
	my $num_in_line = @in_line;

	print "\tConfig line $config_line_num\t".
	    "$num_in_line sections\n" if $verbose;

	# Only try to parse the inline if it has the 
	# expected number of componenets
	if ($num_in_line == 2) {
	    $mask_libs[$i][0] = $in_line[0];
	    $mask_libs[$i][1] = $in_line[1];
	    
	    # Only print for debug runs
	    if ($debug) {
		print STDERR "INLINE SPLIT, i:$i \n";
		print STDERR "\t\t".$in_line[0]."\n";
		print STDERR "\t\t".$in_line[1]."\n";

		print "VAL1 IS:";
		print $mask_libs[$i][0]."\n";
		
		print "VAL2 IS:";
		print $mask_libs[$i][1]."\n";
	    } # End of print for debug runs

	    $i++;
	}
    } # End of unless comment line
} # End of while CONFIGFILE
close CONFIGFILE;

my $num_libs = @mask_libs;
my $num_files = @fasta_files;
my $num_proc_total = $num_libs * $num_files;

#-----------------------------+
# SHOW ERROR IF NO LIBS IN    |
# THE CONFIG FILE             |
#-----------------------------+ 
if ($num_libs == 0) {
    print "\a";
    print STDERR "\nERROR: No library files were indicated in the ".
	"config file\n$config_file\n";
    exit;
}

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

if ($verbose) {
    print STDERR "\n";
    print STDERR "NUM FILES: $num_files\n";
    print STDERR "NUM LIBS:  $num_libs\n";
    print STDERR "NUM PROC:  $num_proc_total\n";
    print STDERR "\n";
}

#-----------------------------+
# SHOW ERROR IF ONE OF THE    |
# MASK LIBS DOES NOT EXIST    |
#-----------------------------+
print STDERR "Checking mask libs ...\n" if $verbose;

for ($i=0; $i<$num_libs; $i++) {

    print "i $i\n";

    $rep_db_name = $mask_libs[$i][0];
    $rep_db_path = $mask_libs[$i][1];

    unless (-e $rep_db_path) {
	print "\a";
	print "\nERROR: The following masking library could not be found:\n".
	    $rep_db_path."\n";
	exit;
    }
}


#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
print "Creating output dir ...\n" unless $quiet;
unless (-e $outdir) {
    mkdir $outdir ||
	die "\aERROR: Could not create the output directory:\n$outdir\n";
}


#-----------------------------+
# RUN REPEAT MASKER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# rep_libs ARRAY              |
#-----------------------------+

for my $ind_file (@fasta_files)
{
    
    $file_num++;

    # Reset search name to null
    $search_name = "";
    
    if ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }
	
    # The following names are the names as produced by
    # RepeatMasker
    $file_to_mask = $indir.$ind_file;
    $repmask_outfile = $file_to_mask.".out";
    $repmask_catfile = $file_to_mask.".cat";
    $repmask_tbl_file = $file_to_mask.".tbl";
    $repmask_masked_file = $file_to_mask.".masked";
    
    $proc_num++;
    print LOG "\n\nProcess $proc_num of $num_proc_total.\n" if $logfile;
    print "\n\n+-----------------------------------------------------------+\n"
	unless $quiet;
    print "| Process $proc_num of $num_proc_total.\n" unless $quiet;
    print "+-----------------------------------------------------------+\n"
	unless $quiet;

    #-----------------------------+
    # SHOW BASIC RUN INFO
    #-----------------------------+
    print "\n";
    print "\tINFILE: $ind_file\n";
    print "\t  ROOT: $name_root\n";

    #-----------------------------+
    # MAKE OUTPUT DIR             |
    #-----------------------------+
    # The base output dir for the BAC
    $bac_out_dir = $outdir.$name_root."/";
    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 

    #-----------------------------+
    # MAKE RM OUTPUT DIR          |
    #-----------------------------+
    # dir for the repeat maske output
    $bac_rep_out_dir = "$bac_out_dir"."rm/";
    mkdir $bac_rep_out_dir, 0777 unless (-e $bac_rep_out_dir); 
    
    #-----------------------------+
    # MAKE GFF OUTPUT DIR         |
    #-----------------------------+
    $gff_out_dir = "$bac_out_dir"."gff/";
    mkdir $gff_out_dir, 0777 unless (-e $gff_out_dir); 

    #-----------------------------+
    # FOR EACH DB IN THE          |
    # mask_libs ARRAY             |
    #-----------------------------+
    for ($i=0; $i<$num_libs; $i++) {
 	
	$rep_db_name = $mask_libs[$i][0];
	$rep_db_path = $mask_libs[$i][1];

	#-----------------------------+
	# GET THE STRING TO SEARCH    |
	# FOR IN THE RM OUTPUT        |
	#-----------------------------+
	# The name used by Repeat Masker is taken from the FASTA header
	# Only the first twenty characters of the FASTA header are used
	open (IN, $file_to_mask);
	while (<IN>) {
	    chomp;
	    if (/^\>(.*)/) {
		print "\tFASTA HEADER:\n\t$_\n" if  $verbose;
		print "\tSEQ_ID:\n\t$1\n" if $verbose;
		
		$search_name = $1;
		
		my $test_len = 20;
		my $cur_len = length ($search_name);
		if ( $cur_len > $test_len ) {
		    $search_name = substr ($_, 0, 20);
		} 
	    } # End of in fasta header file
	}
	close IN;

	#$gff_el_out = $indir.$rep_db_name."_".$ind_file."_EL.gff";
	#$xml_el_out = $indir.$rep_db_name."_".$ind_file."_EL.game.xml"; 
 	#$gff_alldb_out = $indir."ALLDB_".$ind_file.".gff";
	#$xml_alldb_out = $indir."ALLDB_".$ind_file."game.xml";

	# Renamed 09/11/2007
	$gff_el_out = $indir.$name_root."_".$rep_db_name.".gff";
	$xml_el_out = $indir.$name_root."_".$rep_db_name.".game.xml"; 
 	$gff_alldb_out = $indir.$name_root."_ALLDB.gff";
	$xml_alldb_out = $indir.$name_root."_ALLDB.game.xml";
	
	if ($rm_path) {
	    $cmd_repmask = $rm_path.
		" -lib ".$rep_db_path.
		" -pa ".$num_proc.
		" -engine ".$engine.
		" -xsmall".
		" $file_to_mask";
	}
	else {
	    $cmd_repmask = "RepeatMasker".
		" -lib ".$rep_db_path.
		" -pa ".$num_proc.
		" -engine ".$engine.
		" -xsmall".
		" $file_to_mask";
	}       

	#-----------------------------+
	# SHOW THE USER THE COMMANDS  | 
	# THAT WILL BE USED           |
	#-----------------------------+
	if ($verbose) {
	    print "\n";
	    print "+-----------------------------+\n";
	    print "| CONVERT COMMANDS            |\n";
	    print "+-----------------------------+\n";
	    print "\tSEARCH:   ".$search_name."\n";
	    print "\tOUTFILE:  ".$repmask_outfile."\n";
	    print "\tDB-NAME:  ".$rep_db_name."\n";
	    print "\tGFF-FILE: ".$gff_alldb_out."\n";
	    
	    
	    print "\n";
	    print "+-----------------------------+\n";
	    print "| REPEATMASKER COMMANDS       |\n";
	    print "+-----------------------------+\n";
	    print "\tLIB-NAME: ".$rep_db_name."\n";
	    print "\tLIB-PATH: ".$rep_db_path."\n";
	    print "\tEL-OUT:   ".$gff_el_out."\n";
	    print "\tREPCMD:   ".$cmd_repmask."\n";
	    print "\n\n";
	}

	#-----------------------------+
	# PRINT INFO TO LOG FILE      | 
	#-----------------------------+
	if ($logfile) {
	    print LOG "\tLib Name: ".$rep_db_name."\n";
	    print LOG "\tLib Path: ".$rep_db_path."\n";
	    print LOG "\tEL Out:   ".$gff_el_out."\n";
	    print LOG "\tRepCmd:   ".$cmd_repmask."\n";
	    print LOG "\n\n";
	}

	unless ( $test ) {
	    $msg = "\nERROR:\n".
		"Could not complete system cmd\n$cmd_repmask\n";
	    system ( $cmd_repmask );
	}


	unless ( $test ) {
	    rmout_to_gff($repmask_outfile, $gff_el_out, ">");
	}


	unless ( $test ) {
	    rmout_to_gff( $repmask_outfile, $gff_alldb_out, ">>");
	}


	#-----------------------------+
	# CONVERT THE FILES FROM GFF  | 
	# FORMAT TO A MORE USABLE     |
	# FORMAT FOR APOLLO           |
	#-----------------------------+
	if ($apollo) {
	    print "\n\n\n\nCONVERTING GFF TO GAME\n\n\n" if $verbose;
	    &apollo_convert ( $gff_el_out, "gff", $xml_el_out, "game", 
			      $file_to_mask, "none" );  
	}

	#-----------------------------+
	# FILES MOVED TO HERE         |
	#-----------------------------+
	$repmask_out_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".rm.out";
	$repmask_cat_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".rm.cat";
	$repmask_tbl_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".rm.tbl";
	$repmask_masked_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".masked.fasta";
	
	$repmask_el_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".gff";
	$repmask_xml_el_cp = $bac_rep_out_dir.$name_root."_".$rep_db_name.
	    ".game.xml"; 

	#///////////////////////////////////////
	# This is another copy of the masked file
	# this will allow me to put all of the masked files in 
	# a single location and have a shorter name
	# These will all be placed in the $outdir
	$repmask_local_cp = $outdir.$name_root.".masked.fasta";
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



	# THE FOLLOWING ADDED 09/28/2006
	# REmoved 09/11/2007
#	$repmask_log = $indir.$rep_db_name."_".$ind_file.".log";
#	$repmask_log_cp = $bac_rep_out_dir.$rep_db_name."_".$ind_file.".log";
#	$msg = "\nERRORCan not move\n\t".$repmask_log."\n\t".
#	    $repmask_log_cp."\n";
#	move ( $repmask_log, $repmask_log_cp) ||
#	    print STDERR $msg;

	#-----------------------------+
	# MAKE A COPY OF THE MASKED   |
	# FASTA FILE TO A SINGLE DIR  |
	#-----------------------------+
	$msg = "\nERROR: Can not copy ".$repmask_masked_file." to\n".
	    $repmask_local_cp."\n";
	copy ( $repmask_masked_file, $repmask_local_cp  ) ||
	    print STDERR $msg;

	#-----------------------------+
	# MOVE THE RM OUTPUT FILES TO |
	# THE TARGET DIR	      |
	#-----------------------------+
	$msg = "\nERROR: Can not move ".$repmask_outfile.
	    " to\n ".$repmask_cat_cp."\n";
	move ( $repmask_outfile, $repmask_cat_cp ) ||
	    print STDERR $msg;

	$msg = "\nERROR: Can not move ".$repmask_catfile."\n";
	move ( $repmask_catfile, $repmask_cat_cp ) ||
	    print STDERR $msg;
	
#	$msg = "\nERROR: The table file could not be moved from".
#	    "$repmask_tbl_file to $repmask_tbl_cp";
#	move ( $repmask_tbl_file, $repmask_tbl_cp ) ||
#	    print STDERR $msg;    
	
	$msg = "\nERROR: Can not move ".$repmask_masked_file."\n";
	move ( $repmask_masked_file, $repmask_masked_cp ) ||
	    print STDERR $msg;

	$msg = "\nERROR: Can not move ".$gff_el_out."\n";
	move ( $gff_el_out , $repmask_el_cp ) ||
	    print STDERR $msg;
	
	if ($apollo) {
	    $msg = "\nERROR: Can not move ".$xml_el_out."\n";
	    move ( $xml_el_out, $repmask_xml_el_cp ) ||
		print STDERR $msg;
	}

    } # End of for LibData
    
    #-----------------------------+ 
    # THE APOLLO CONVERT FOR THE  |
    # ENTIRE SET OF REPEAT LIBS   |
    # FOR A GIVEN SEQUENCE FILE   |
    # THAT IS BEING MASKED.       |
    #-----------------------------+
    if ($apollo) {
	apollo_convert ( $gff_alldb_out, "gff", $xml_alldb_out , "game", 
			  $file_to_mask, "none" );  
    }
    
    my $repmask_all_gff_cp = $bac_rep_out_dir.$name_root."_ALLDB.rm.gff";
    my $repmask_all_xml_cp = $bac_rep_out_dir.$name_root."_ALLDB.rm.game.xml";

    $msg = "\nCan not move ".$gff_alldb_out."\n";
    move ( $gff_alldb_out, $repmask_all_gff_cp ) ||
	print STDERR $msg;
    
    if ($apollo) {
	$msg = "\nCan not move ".$xml_alldb_out."\n";
	move ( $xml_alldb_out, $repmask_all_xml_cp ) ||
	    print STDERR $msg;
    }


    # TEMP EXIT FOR DEBUG, WIll JUST RUN FIRST FILE TO BE MASKED
    if ($debug) {
	if ($file_num > $file_num_max ) {
	    print "\nDebug run finished\n\n";
	    exit;
	}
    }
    
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
	print RM_OUT "$rmout[4]\t";            # qry sequence name
	print RM_OUT "repeatmasker:trep9\t";   # software used
	print RM_OUT "exon\t";                 # attribute name
	print RM_OUT "$rmout[5]\t";            # start
	print RM_OUT "$rmout[6]\t";            # stop
	print RM_OUT "$rmout[0]\t";            # smith waterman score"
	print RM_OUT "+\t";                    # Postive strand
	print RM_OUT ".\t";                    # frame
	print RM_OUT "$rmout[9]";              # attribute
	print RM_OUT "\n";

	#print RM_OUT "$rmout[4]\t";            # qry sequence name
	#print RM_OUT "repeatmasker:trep9\t";   # software used
	#print RM_OUT "exon\t";                 # attribute name
	#print RM_OUT "$rmout[5]\t";            # start
	#print RM_OUT "$rmout[6]\t";            # stop
	#print RM_OUT "$rmout[0]\t";            # smith waterman score"
	#print RM_OUT "-\t";                    # Negative strand
	#print RM_OUT ".\t";                    # frame
	#print RM_OUT "$rmout[9]";              # attribute
	#print RM_OUT "\n";

    }

    return;
    
}

=head1 NAME

batch_mask.pl - Run RepeatMasker and parse results to a gff format file. 

=head1 VERSION

This documentation refers to program version $Rev$

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

VERSION: $Rev$

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
