#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_cnv_ta2ap.pl - Convert TriAnnotation to Apollo GFF  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/10/2006                                       |
# UPDATED: 07/18/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the TriAnnotation gff output to a format that   |
#  can be used by the Apollo Genome Annotation Curation     |
#  program.                                                 |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

batch_cnv_ta2ap.pl - Convert TriAnnotation to Apollo GFF

=head1 VERSION

This documentation refers to batch_mask version 1.0

=head1 SYNOPSIS

 Usage:
 batch_mask.pl -i DirToProcess -o OutDir

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

=item -p,--num-proc

The number of processors to use for RepeatMasker. Default is one.

=item --engine

The repeatmasker engine to use: [crossmatch|wublast|decypher].
The default is to use crossmatch.

=item --apollo

Use the apollo program to convert the file from gff to game xml.
The default is not to use apollo.

=item --rm-path

The full path to the RepeatMasker binary.

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

File::Copy

=item *

Getopt::Long

=back

=head1 BUGS AND LIMITATIONS

=head2 TO DO

=over 2

=item *

Load the RepLibs array from a config file.

=item *

Make the results compatable for an upload to a chado
database.

=item *

Make it a variable to possible to put the gff output (1) all in positive 
strand, (2) all in negative strand, (3) alignment to positive or
negative strand, (4) cumulative in both positive and negative strand.
Current behavior will be to do number 4 above.

=back

=head2 Limitations

=over

=item *

Currently must use short names in the FASTA file.

=item *

This program has been tested with RepeatMasker v  3.1.6

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
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file

my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc
my $rm_path;                   # Full path to the repeatmasker binary
my $ap_path;                   # Full path to the apollo program

# Vars with default values
my $engine = "crossmatch";

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
my $num_proc = 1;              # Number of processors to use

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "rm-path=s"    => \$rm_path,
		    "ap-path=s",   => \$ap_path,
		    "logfile=s"    => \$logfile,
		    "p|num-proc=s" => \$num_proc,
		    "engine=s"     => \$engine,
		    # Booleans
		    "apollo"       => \$apollo,
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);

#-----------------------------+
# HARD CODE VARIABLES         |
#-----------------------------+
# TempWorkDir - A directory to temporarily store the data that will later
#               be transfered to a more appropriate long term storage
#               location. This could be a place to put the original fasta file
#               and associated intermediate output.
# FileToMask - The fasta file that will be masked using RepeatMasker
# GFFOutfile - The path to write the gff outfile
# RepMaskOutfile - The path of the *.out file produced by RepeatMasker
# SeqName - The name of the sequence file in the original [name].fasta
#           file that was used as input into RepeatMasker
# RepDb - The short name of the repeat mask DB that was used

# Many of the following variables will need to be changed to set to 
# array elements, but these are hard coded for now to allow for testing and
# getting the program up and running.

# MakeGffDBCmd - Command to make a GFF file with the data grouped by the 
#                the repeat database that was used for assignment.
# MakeGffElCmd - Command to make a GFF file with the data grouped by the
#                type of element that was identified.



#my $bac_parent_dir = "/scratch/jestill/wheat/";  

#my $local_cp_dir = "/scratch/jestill/wheat/copy/"; 
#my $LocalCpDir = $local_cp_dir;

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
# This is an array of arrays:
# The nested array should contain the following info:
# [0] - Name of the repeat library
#       This will be used to name the output files from the analysis
#       and to name the data tracks that will be used by Apollo.
# [1] - Path of the repeat library
# The number of records in the fasta file shown in brackets
# afte the description of the db
@mask_libs = (
    #-----------------------------+
    # TREP v 9                    |
    #-----------------------------+
    ["TREP9",  
     "/db/jlblab/repeats/TREP9.nr.fasta"]
    );

my $num_libs = @mask_libs;
my $num_files = @fasta_files;
my $num_proc_total = $num_libs * $num_files;

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
# SHOW ERROR IF ONE OF THE    |
# MASK LIBS DOES NOT EXIST    |
#-----------------------------+
print "Checking mask libs ...\n" unless $quiet;
for $ind_lib (@mask_libs) {

    $RepDbName = @$ind_lib[0];
    $RepDbPath = @$ind_lib[1];

    unless (-e $RepDbPath) {
	print "\a";
	print "\nERROR: The following masking library could not be found:\n".
	    $RepDbPath."\n";
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
	
    $file_to_mask = $indir.$ind_file;
    $RepMaskOutfile = $file_to_mask.".out";
    $RepMaskCatFile = $file_to_mask.".cat";
    $RepMaskTblFile = $file_to_mask.".tbl";
    $RepMaskMaskedFile = $file_to_mask.".masked";
    
    $ProcNum++;
    print LOG "\n\nProcess $ProcNum of $num_proc_total.\n" if $logfile;
    print "\n\n+-----------------------------------------------------------+\n"
	unless $quiet;
    print "| Process $ProcNum of $num_proc_total.\n" unless $quiet;
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
    $bac_out_dir = $outdir.$name_root."/";
    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 

    for $ind_lib (@mask_libs)
    {
	
	$RepDbName = @$ind_lib[0];
	$RepDbPath = @$ind_lib[1];

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
		print "\tINSIDE:\n\t$1\n" if $verbose;
		
		$search_name = $1;
		
		my $test_len = 20;
		my $cur_len = length ($search_name);
		if ( $cur_len > $test_len ) {
		    $search_name = substr ($_, 0, 20);
		} 
	    } # End of in fasta header file
	}
	close IN;

	$GffElOut = $indir.$RepDbName."_".$ind_file."_EL.gff";
	$XmlElOut = $indir.$RepDbName."_".$ind_file."_EL.game.xml"; 
 	$GffAllDbOut = $indir."ALLDB_".$ind_file.".gff";
	$XmlAllDbOut = $indir."ALLDB_".$ind_file."game.xml";
	

	if ($rm_path) {
	    $RepMaskCmd = $rm_path.
		" -lib ".$RepDbPath.
		" -pa ".$num_proc.
		" -engine ".$engine.
		" -xsmall".
		" $file_to_mask";
	}
	else {
	    $RepMaskCmd = "RepeatMasker".
		" -lib ".$RepDbPath.
		" -pa ".$num_proc.
		" -engine ".$engine.
		" -xsmall".
		" $file_to_mask";
	}
       


	#-----------------------------+
	# SHOW THE USER THE COMMANDS  | 
	# THAT WILL BE USED           |
	#-----------------------------+

	print "\n";
	print "+-----------------------------+\n";
	print "| CONVERT COMMANDS            |\n";
	print "+-----------------------------+\n";
	print "\tSEARCH:   ".$search_name."\n";
	print "\tOUTFILE:  ".$RepMaskOutfile."\n";
	print "\tDB-NAME:  ".$RepDbName."\n";
	print "\tGFF-FILE: ".$GffAllDbOut."\n";


	print "\n";
	print "+-----------------------------+\n";
	print "| REPEATMASKER COMMANDS       |\n";
	print "+-----------------------------+\n";
	print "\tLIB-NAME: ".$RepDbName."\n";
	print "\tLIB-PATH: ".$RepDbPath."\n";
	print "\tEL-OUT:   ".$GffElOut."\n";
	print "\tREPCMD:   ".$RepMaskCmd."\n";
	print "\n\n";

	#-----------------------------+
	# PRINT INFO TO LOG FILE      | 
	#-----------------------------+
	if ($logfile) {
	    print LOG "\tLib Name: ".$RepDbName."\n";
	    print LOG "\tLib Path: ".$RepDbPath."\n";
	    print LOG "\tEL Out:   ".$GffElOut."\n";
	    print LOG "\tRepCmd:   ".$RepMaskCmd."\n";
	    print LOG "\n\n";
	}

	# Turned off while working 07/13/2007
	unless ( $test ) {
	    $msg = "\nERROR:\n".
		"Could not complete system cmd\n$RepMaskCmd\n";
	    system ( $RepMaskCmd );
	}


	unless ( $test ) {
	    rmout_to_gff($RepMaskOutfile, $GffElOut, ">");
	}


	unless ( $test ) {
	    rmout_to_gff( $RepMaskOutfile, $GffAllDbOut, ">>");
	}


	#-----------------------------+
	# CONVERT THE FILES FROM GFF  | 
	# FORMAT TO A MORE USABLE     |
	# FORMAT FOR APOLLO           |
	#-----------------------------+

	# APOLLO FUNCTION WILL NOT WORK ON ALTIX
	# SO THIS HAS BEEN COMMENTED OUT
	#print "\n\n\n\nCONVERTING\n\n\n";
	if ($apollo) {
	    &apollo_convert ( $GffElOut, "gff", $XmlElOut, "game", 
			      $file_to_mask, "none" );  
	}

	#-----------------------------+
	# COPY THE RM OUTPUT FILES TO | 
	# THE RM (REPEATMASK) FOLDER  |
	# MAKE THE DIR IF NEEDED      |
	#-----------------------------+
	# dir for the repeat maske output
	$bac_rep_out_dir = "$bac_out_dir"."rm/";
	mkdir $bac_rep_out_dir, 0777 unless (-e $bac_rep_out_dir); 

	#-----------------------------+
	# FILES MOVED TO HERE         |
	#-----------------------------+
	$RepMaskOutCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    ".rm.out";
	$RepMaskCatCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    ".rm.cat";
	$RepMaskTblCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    ".rm.tbl";
	$RepMaskMaskedCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    ".masked.fasta";
	$RepMaskElCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    "_EL.gff";
	$RepMaskXmlElCp = $bac_rep_out_dir.$name_root."_".$RepDbName.
	    "_EL.game.xml"; 

	#///////////////////////////////////////
	# This is another copy of the masked file
	# this will allow me to put all of the masked files in 
	# a single location and have a shorter name
	# These will all be placed in the $outdir
	$RepMaskLocalCp = $outdir.$name_root.".masked.fasta";
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	# THE FOLLOWING ADDED 09/28/2006
	$RepMaskLog = $indir.$RepDbName."_".$ind_file.".log";
	$RepMaskLogCp = $bac_rep_out_dir.$RepDbName."_".$ind_file.".log";
	$msg = " Can not move\n\t".$RepMaskLog."\n\t".
	    $RepMaskLogCp."\n";
	move ( $RepMaskLog, $RepMaskLogCp) ||
	    print LOG $msg if $logfile;

	#-----------------------------+
	# MAKE A COPY OF THE MASKED   |
	# FASTA FILE TO A SINGLE DIR  |
	#-----------------------------+
	$msg = "Can not copy ".$RepMaskMaskedFile." to\n".
	    $RepMaskLocalCp."\n";
	copy ( $RepMaskMaskedFile, $RepMaskLocalCp  ) ||
	    print LOG $msg if $logfile;

	#-----------------------------+
	# MOVE THE RM OUTPUT FILES TO |
	# THE TARGET DIR	      |
	#-----------------------------+
	$msg = "Can not move ".$RepMaskOutfile." to\n ".$RepMaskOutCp."\n";
	move ( $RepMaskOutfile, $RepMaskOutCp ) ||
	    print LOG $msg if $logfile;

	$msg = "Can not move ".$RepMaskCatFile."\n";
	move ( $RepMaskCatFile, $RepMaskCatCp ) ||
	    print LOG $msg if $logfile;
	
	$msg = "The table file could not be moved from".
	    "$RepMaskTblFile to $RepMaskTblCp";
	move ( $RepMaskTblFile, $RepMaskTblCp ) ||
	    print LOG $msg if $logfile;    
	
	$msg = "Can not move ".$RepMaskMaskedFile."\n";
	move ( $RepMaskMaskedFile, $RepMaskMaskedCp ) ||
	    print LOG $msg if $logfile;

	$msg = "Can not move ".$GffElOut."\n";
	move ( $GffElOut , $RepMaskElCp ) ||
	    print LOG $msg if $logfile;
	
	if ($apollo) {
	    $msg = "Can not move ".$XmlElOut."\n";
	    move ( $XmlElOut, $RepMaskXmlElCp ) ||
		print LOG $msg if $logfile;
	}

    } # End of for LibData
    
    #-----------------------------+ 
    # THE APOLLO CONVERT FOR THE  |
    # ENTIRE SET OF REPEAT LIBS   |
    # FOR A GIVEN SEQUENCE FILE   |
    # THAT IS BEING MASKED.       |
    #-----------------------------+
    if ($apollo) {
	apollo_convert ( $GffAllDbOut, "gff", $XmlAllDbOut , "game", 
			  $file_to_mask, "none" );  
    }
    
    $RepMaskALL_GFFCp = $bac_rep_out_dir.$name_root."_ALLDB.rm.gff";
    $RepMaskAll_XMLCp = $bac_rep_out_dir.$name_root."_ALLDB.rm.game.xml";

    $msg = "Can not move ".$GffAllDbOut."\n";
    move ( $GffAllDbOut, $RepMaskALL_GFFCp ) ||
	print LOG $msg if $logfile;
    
    if ($apollo) {
	$msg = "Can not move ".$XmlAllDbOut."\n";
	move ( $XmlAllDbOut, $RepMaskAll_XMLCp ) ||
	    print LOG $msg if $logfile;
    }


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
	"  batch_mask.pl -i DirToProcess -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --rm-path      # Full path to repeatmasker binary".
	"  --engine       # The repeatmasker engine to use:\n".
	"                 # [crossmatch|wublast|decypher]\n".
	"                 # default is to use crossmatch\n".
	"  --num-proc     # Number of processors to use for RepeatMasker\n".
	"                 # default is one.\n".
	"  --apollo       # Convert output to game.xml using apollo\n".
	"                 # default is not to use apollo\n".
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
